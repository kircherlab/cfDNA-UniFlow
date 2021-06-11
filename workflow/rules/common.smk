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

### not fully implemented -> cases: fastQ files as input -> SR,PE 
def get_NGmerge_input(wildcards):
    sample=wildcards.SAMPLE
    inpath=samples.loc[sample].loc["path"]
    if ".bam" in inpath.lower():
        return {"r1":"results/{ID}/fastq/{SAMPLE}_R1.fastq.gz","r2":"results/{ID}/fastq/{SAMPLE}_R2.fastq.gz"}
