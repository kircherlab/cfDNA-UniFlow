from snakemake.utils import validate
import pandas as pd
from pathlib import Path
import sys

##### load config and sample sheets #####


configfile: "config/config.yaml"


# validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
# validate(samples, schema="../schemas/samples.schema.yaml")


def get_final_output():
    final_output = []

    final_output.extend(
        expand(
            "results/{ID}/flagstat/{SAMPLE}_flagstat.txt.gz",
            zip,
            ID=samples["ID"],
            SAMPLE=samples["sample"]
        )
    )
    final_output.extend(
        expand(
            "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
            zip,
            ID=samples["ID"],
            SAMPLE=samples["sample"],
            GENOME=samples["genome_build"]
        )
    )
    final_output.extend(
        expand(
            "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
            zip,
            ID=samples["ID"],
            SAMPLE=samples["sample"],
            GENOME=samples["genome_build"]
        )
    )
    final_output.extend(
        expand(
            "results/{ID}/qc/multiqc.html",
            ID=samples["ID"].unique()
        )
    )

    return final_output

def get_read_group(sample):
    ID = samples.loc[sample].loc["ID"]
    library = samples.loc[sample].loc["library_name"]
    platform = samples.loc[sample].loc["platform"]
    info = samples.loc[sample].loc["info"]
    if pd.isna(info):
        info=""
    else:    
        info="_"+info
    RGID = f"{ID}_{sample}{info}"
    RG = f"@RG\\tID:{sample}\\tSM:{RGID}\\tLB:{library}\\tPL:{platform}"
    return RG

### not fully implemented -> cases: fastQ files as input -> SE,PE 
def get_trimming_input(wildcards):
    sample=wildcards.SAMPLE
    inpath=samples.loc[sample].loc["path"]
    if ".bam" in inpath.lower():
        return {"r1":"results/{ID}/fastq/{SAMPLE}_R1.fastq.gz","r2":"results/{ID}/fastq/{SAMPLE}_R2.fastq.gz"}


def get_reference(wildcards):
    genome_build = wildcards.GENOME
    gpath = config[genome_build]["reference"]
    p=Path(gpath)
    if p.exists() and not p.is_dir():
        try:
            p.open().close()
            return p.as_posix()
        except PermissionError as f:
            print(f,file = sys.stderr)
            print(f"Please check file permissions of {gpath} or remove path from config.yaml \
            download a reference.",file = sys.stderr)
            raise f
    else:
        return f"resources/reference/{genome_build}.fa"

def get_ref_url(wildcards):
    genome_build = wildcards.GENOME
    url_dict={
        "hg19":"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
        "hg38":"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    }
    return url_dict[genome_build]

