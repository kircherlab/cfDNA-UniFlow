from snakemake.utils import validate
import pandas as pd
from pathlib import Path
import sys

##### load config and sample sheets #####


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

regions = pd.read_csv(config["regions"], sep="\t").set_index("target", drop=False)
regions.index.names = ["region_id"]
validate(regions, schema="../schemas/regions.schema.yaml")


def get_final_output():
    final_output = []

    final_output.extend(
        expand(
            "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
            zip,
            ID=samples["ID"],
            SAMPLE=samples["sample"],
            GENOME=samples["genome_build"],
        )
    )
    final_output.extend(
        expand(
            "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
            zip,
            ID=samples["ID"],
            SAMPLE=samples["sample"],
            GENOME=samples["genome_build"],
        )
    )

    final_output.extend(
        expand("results/{ID}/qc/multiqc.html", ID=samples["ID"].unique())
    )

    if config["utility"]["GCbias-correction"]:
        final_output.extend(
            expand(
                expand(
                    "results/{ID}/signals/signal-corrected/{target_region}.{SAMPLE}-corrected_WPS.{GENOME}.csv.gz",
                    zip,
                    ID=samples["ID"],
                    SAMPLE=samples["sample"],
                    GENOME=samples["genome_build"],
                    allow_missing=True,
                ),
                target_region=regions["target"],
            )
        )

    if config["utility"]["GCbias-plot"]:
        final_output.extend(
            expand(
                expand(
                    "results/{ID}/GCBias/plots/{SAMPLE}-GCbias-plot_{blacklist}.{GENOME}.png",
                    zip,
                    ID=samples["ID"],
                    SAMPLE=samples["sample"],
                    GENOME=samples["genome_build"],
                    allow_missing=True,
                ),
                blacklist=["repeatmasker"],
            )
        )

        final_output.extend(
            set(
                expand(
                    expand(
                        "results/{ID}/signals/GCcorrection-plots/{target_region}.{status_name}-GCcorrected_{signal}.{GENOME}.png",
                        zip,
                        ID=samples["ID"],
                        GENOME=samples["genome_build"],
                        status_name=samples["status"],
                        allow_missing=True,
                    ),
                    signal="COV" if config["signal"].lower() == "coverage" else "WPS",
                    target_region=regions["target"],
                )
            )
        )

    if config["utility"]["ichorCNA"]:
        final_output.extend(
            expand(
                "results/{ID}/icorCNA/readcounts/{SAMPLE}_processed.{GENOME}.wig",
                zip,
                ID=samples["ID"],
                SAMPLE=samples["sample"],
                GENOME=samples["genome_build"],
            )
        )
        final_output.extend(
            expand(
                "results/{ID}/icorCNA/{SAMPLE}_processed_{GENOME}",
                zip,
                ID=samples["ID"],
                SAMPLE=samples["sample"],
                GENOME=samples["genome_build"],
            )
        )

    if config["utility"]["case-control-plot"]:
        final_output.extend(
            set(
                expand(
                    expand(
                        "results/{ID}/signals/case-control-plots/{target_region}.{case_name}-vs-{control_name}-{correction}_{signal}.{GENOME}.png",
                        zip,
                        ID=samples["ID"],
                        GENOME=samples["genome_build"],
                        case_name=samples["status"].loc[
                            samples["status"] != config["control_name"]
                        ],
                        allow_missing=True,
                    ),
                    control_name=config["control_name"],
                    signal="COV" if config["signal"].lower() == "coverage" else "WPS",
                    correction="corrected"
                    if config["utility"]["GCbias-correction"]
                    else "uncorrected",
                    target_region=regions["target"],
                )
            )
        )

    return final_output


def get_read_group(sample):
    ID = samples.loc[sample].loc["ID"]
    library = samples.loc[sample].loc["library_name"]
    platform = samples.loc[sample].loc["platform"]
    info = samples.loc[sample].loc["info"]
    if pd.isna(info):
        info = ""
    else:
        info = "_" + info
    RGID = f"{ID}_{sample}{info}"
    RG = f"@RG\\tID:{sample}\\tSM:{RGID}\\tLB:{library}\\tPL:{platform}"
    return RG


### get input files for trimming rules (NGmerge/trimmomatic) based on input format
def get_trimming_input(wildcards):
    sample = wildcards.SAMPLE
    bam = samples.loc[sample].loc["bam"]
    fq1 = samples.loc[sample].loc["fq1"]
    fq2 = samples.loc[sample].loc["fq2"]

    if ".fastq" in fq1.lower() and ".fastq" in fq2.lower():
        return {
            "r1": fq1,
            "r2": fq2,
        }
    elif ".bam" in bam.lower():
        return {
            "r1": "results/{ID}/fastq/{SAMPLE}_R1.fastq.gz",
            "r2": "results/{ID}/fastq/{SAMPLE}_R2.fastq.gz",
        }
    else:
        raise ValueError(
            f"No raw fastq files or bam file found for sample: {sample} (bam: {bam}; fq1: {fq1}; fq2: {fq2})."
            f"Please check the your sample sheet."
        )


def get_SEtrimming_input(wildcards):
    sample = wildcards.SAMPLE
    bam = samples.loc[sample].loc["bam"]
    fq1 = samples.loc[sample].loc["fq1"]
    fq2 = samples.loc[sample].loc["fq2"]

    if ".bam" in bam.lower():
        return ["results/{ID}/fastq/{SAMPLE}_single_read.fastq.gz"]
    elif ".fastq" in fq1.lower() and ".fastq" in fq2.lower():
        return ["results/{ID}/fastq/{SAMPLE}_single_read.fastq.gz"]
    elif ".fastq" in fq1.lower() and not ".fastq" in fq2.lower():
        return [fq1]
    elif ".fastq" in fq2.lower() and not ".fastq" in fq1.lower():
        return [fq2]
    else:
        raise ValueError(
            f"Something went wrong for sample: {sample} (bam: {bam}; fq1: {fq1}; fq2: {fq2})."
            f"Please check the your sample sheet."
        )


### get reference based on genome_build provided in samples.tsv
def get_reference(wildcards):
    genome_build = wildcards.GENOME
    gpath = config[genome_build]["reference"]
    p = Path(gpath)
    if p.exists() and not p.is_dir():
        try:
            p.open().close()
            return p.as_posix()
        except PermissionError as f:
            print(
                f"Please check file permissions of {gpath} or remove path from config.yaml \
            download a reference.",
                file=sys.stderr,
            )
            raise f
    else:
        return f"resources/reference/{genome_build}.fa"


### provides download URLs for UCSC human reference files in fasta format
def get_fasta_ref_url(wildcards):
    genome_build = wildcards.GENOME
    url_dict = {
        "hg19": "rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz",
        "hg38": "rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz",
    }
    return url_dict[genome_build]


### provides download URLs for UCSC human reference files in twobit format
def get_twobit_ref_url(wildcards):
    genome_build = wildcards.GENOME
    url_dict = {
        "hg19": "rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.2bit",
        "hg38": "rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.2bit",
    }
    return url_dict[genome_build]


### returns properly formattet trimming rules based on config.yaml
def get_trimmomatic_trimmers():
    default_adapter = {
        "nexterape-pe.fa": "resources/adapter/NexteraPE-PE.fa",
        "truseq2-pe.fa": "resources/adapter/TruSeq2-PE.fa",
        "truseq2-se.fa": "resources/adapter/TruSeq2-SE.fa",
        "truseq3-pe-2.fa": "resources/adapter/TruSeq3-PE-2.fa",
        "truseq3-pe.fa": "resources/adapter/TruSeq3-PE.fa",
        "truseq3-se.fa": "resources/adapter/TruSeq3-SE.fa",
    }
    trimmers = list()
    t_conf = config["trimmers"]

    if default_adapter.get(t_conf["Illuminaclip"]["adapter_files"].lower()) is not None:
        adapter_path = Path(
            default_adapter.get(t_conf["Illuminaclip"]["adapter_files"].lower())
        )
    else:
        adapter_path = Path(t_conf["Illuminaclip"]["adapter_files"])

    if isinstance(t_conf["Illuminaclip"]["seedMismatches"], int):
        mismatches = abs(t_conf["Illuminaclip"]["seedMismatches"])
    else:
        mismatches = 4
        print(
            f"Invalid config, setting default value for seedMismatches: {mismatches}",
            file=sys.stderr,
        )

    if isinstance(t_conf["Illuminaclip"]["palindromeClipThreshold"], int):
        read_overlap = abs(t_conf["Illuminaclip"]["palindromeClipThreshold"])
    else:
        read_overlap = 30
        print(
            f"Invalid config, setting default value for palindromeClipThreshold: {read_overlap}",
            file=sys.stderr,
        )

    if isinstance(t_conf["Illuminaclip"]["simpleClipThreshold"], int):
        minimal_overlap = abs(t_conf["Illuminaclip"]["simpleClipThreshold"])
    else:
        minimal_overlap = 10
        print(
            f"Invalid config, setting default value for simpleClipThreshold: {minimal_overlap}",
            file=sys.stderr,
        )

    if isinstance(t_conf["Illuminaclip"]["minAdapterLength"], int):
        minAdapterLength = abs(t_conf["Illuminaclip"]["minAdapterLength"])
    else:
        minAdapterLength = 8
        print(
            f"Invalid config, setting default value for minAdapterLength: {minAdapterLength}",
            file=sys.stderr,
        )
    if t_conf["Illuminaclip"]["keepBothReads"] is True:
        keepBothReads = "True"
    else:
        keepBothReads = "False"

    if adapter_path.exists() and not adapter_path.is_dir():
        try:
            adapter_path.open().close()
            illumina_string = f"ILLUMINACLIP:{adapter_path}:{mismatches}:{read_overlap}:{minimal_overlap}:{minAdapterLength}:{keepBothReads}"
            trimmers.append(illumina_string)
        except PermissionError as f:
            print(
                f"Please check file permissions of {adapter_path} or select one of the provided files.",
                file=sys.stderr,
            )
            raise f

    if (
        t_conf["SLIDINGWINDOW"]["window"] > 0
        and t_conf["SLIDINGWINDOW"]["quality_threshold"] > 0
    ):
        window = t_conf["SLIDINGWINDOW"]["window"]
        threshold = t_conf["SLIDINGWINDOW"]["quality_threshold"]
        trimmers.append(f"SLIDINGWINDOW:{window}:{threshold}")

    if t_conf["LEADING"] > 0:
        LEADING = t_conf["LEADING"]
        trimmers.append(f"LEADING:{LEADING}")
    if t_conf["TRAILING"] > 0:
        TRAILING = t_conf["TRAILING"]
        trimmers.append(f"TRAILING:{TRAILING}")
    if config["length-filter"]["MINLEN"] > 0:
        MINLEN = t_conf["MINLEN"]
        trimmers.append(f"MINLEN:{MINLEN}")

    return trimmers


def get_trimmomatic_extra():
    """Checks if phred-qualty-encoding is set to phred-33 or phred-64."""
    trimmomatic_extra = ""
    if config["phred-quality-encoding"] == "phred-33":
        trimmomatic_extra = "-phred33"
    elif config["phred-quality-encoding"] == "phred-64":
        trimmomatic_extra = "-phred64"

    return trimmomatic_extra


def get_mapping_input(wildcards):
    sample = wildcards.SAMPLE
    bam = samples.loc[sample].loc["bam"]
    fq1 = samples.loc[sample].loc["fq1"]
    fq2 = samples.loc[sample].loc["fq2"]
    mapping_input = dict()
    # Paired End parameters
    trimming_algorithm = config["PE_trimming_algorithm"]
    unmerged = config["mapping"]["paired_end"]["unmerged"]
    singleton = config["mapping"]["paired_end"]["singleton"]
    # Single End parameter
    SEreads = config["mapping"]["single_end"]["SEreads"]

    if trimming_algorithm.lower() == "ngmerge":
        # these options are for Paired End sequencing libraries
        if ".bam" in bam.lower() or (
            ".fastq" in fq1.lower() and ".fastq" in fq2.lower()
        ):
            mapping_input[
                "reads"
            ] = "results/{ID}/NGmerge/merged/{SAMPLE}_merged.filtered.fastq.gz"
            if unmerged:
                mapping_input[
                    "noadapter_R1"
                ] = "results/{ID}/NGmerge/nonmerged/{SAMPLE}_noadapters_1.filtered.fastq.gz"
                mapping_input[
                    "noadapter_R2"
                ] = "results/{ID}/NGmerge/nonmerged/{SAMPLE}_noadapters_2.filtered.fastq.gz"
            if singleton and ".bam" in bam.lower():
                mapping_input[
                    "PEsingleton"
                ] = "results/{ID}/fastq/{SAMPLE}_PEsingleton.filtered.fastq.gz"
        # this option is for reads in Single End sequencing libraries
        if SEreads and (".fastq" in fq1.lower() and not ".fastq" in fq2.lower()):
            mapping_input[
                "reads"
            ] = "results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.filtered.fastq.gz"
        elif SEreads and (".fastq" in fq2.lower() and not ".fastq" in fq1.lower()):
            mapping_input[
                "reads"
            ] = "results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.filtered.fastq.gz"
        elif SEreads and not (".fastq" in fq1.lower() and ".fastq" in fq2.lower()):
            mapping_input[
                "single_reads"
            ] = "results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.filtered.fastq.gz"

    elif trimming_algorithm.lower() == "trimmomatic":
        # these options are for Paired End sequencing libraries
        if ".bam" in bam.lower() or (
            ".fastq" in fq1.lower() and ".fastq" in fq2.lower()
        ):
            mapping_input["reads"] = [
                "results/{ID}/trimmed/trimmomatic/{SAMPLE}.1.fastq.gz",
                "results/{ID}/trimmed/trimmomatic/{SAMPLE}.2.fastq.gz",
            ]
            if unmerged:
                mapping_input[
                    "noadapter_R1"
                ] = "results/{ID}/trimmed/trimmomatic/{SAMPLE}.1.unpaired.fastq.gz"
                mapping_input[
                    "noadapter_R2"
                ] = "results/{ID}/trimmed/trimmomatic/{SAMPLE}.2.unpaired.fastq.gz"
        # this option is for reads in Single End sequencing libraries
        if SEreads and (".fastq" in fq1.lower() and not ".fastq" in fq2.lower()):
            mapping_input[
                "reads"
            ] = "results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.filtered.fastq.gz"
        elif SEreads and (".fastq" in fq2.lower() and not ".fastq" in fq1.lower()):
            mapping_input[
                "reads"
            ] = "results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.filtered.fastq.gz"
        elif SEreads and not (".fastq" in fq1.lower() and ".fastq" in fq2.lower()):
            mapping_input[
                "single_reads"
            ] = "results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.filtered.fastq.gz"
    return mapping_input


def get_mark_duplicates_input(wildcards):
    mark_duplicates_input = dict()
    GCcorrection = config["utility"]["GCbias-correction"]
    blacklist = "repeatmasker"
    if GCcorrection:
        mark_duplicates_input[
            "mapped_reads"
        ] = "results/{{ID}}/corrected_reads/{{SAMPLE}}_GCcorrected_{blacklist}.{{GENOME}}.bam".format(
            blacklist=blacklist,
        )
    else:
        mark_duplicates_input[
            "mapped_reads"
        ] = "results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam"

    return mark_duplicates_input
