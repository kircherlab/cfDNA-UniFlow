from snakemake.utils import validate
import pandas as pd

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
            "results/fastq/{SAMPLE}_{paired_end}.fastq.gz",
            SAMPLE=samples["sample"],
            paired_end=["R1", "R2"],
        )
    )
    final_output.extend(
        expand(
            "results/NGmerge/{SAMPLE}_merged.fastq.gz",
            SAMPLE=samples["sample"]
        )
    )
    final_output.extend(
        expand(
            "results/mapped_reads/{SAMPLE}_all.{GENOME}.bam",
            zip,
            SAMPLE=samples["sample"],
            GENOME=samples["genome_build"]
        )
    )
    final_output.extend(
        expand(
            "results/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
            zip,
            SAMPLE=samples["sample"],
            GENOME=samples["genome_build"]
        )
    )

    return final_output

def get_read_group(sample):
    ID = samples.loc[sample].loc["ID"]
    library = samples.loc[sample].loc["library_name"]
    platform = samples.loc[sample].loc["platform"]
    RGID = f"{ID}_{sample}"
    RG = f"@RG\\tID:{RGID}\\tSM:{sample}\\tLB:{library}\\tPL:{platform}"
    return RG


def get_map_reads_input():
    pass
