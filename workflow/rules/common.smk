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

    return final_output


def get_map_reads_input():
    pass
