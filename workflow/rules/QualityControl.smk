
rule fastqc:
    input:
        "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
    output:
        html="results/{ID}/qc/fastqc/{SAMPLE}_processed.{GENOME}.html",
        zip="results/{ID}/qc/fastqc/{SAMPLE}_processed.{GENOME}_fastqc.zip",
    log:
        "results/logs/{ID}/fastqc/{SAMPLE}_all.{GENOME}.log",
    wrapper:
        "0.75.0/bio/fastqc"


rule samtools_stats:
    input:
        "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
    output:
        "results/{ID}/qc/samtools-stats/{SAMPLE}_processed.{GENOME}.txt",
    log:
        "results/logs/{ID}/fastqc/{SAMPLE}.{GENOME}.log",
    wrapper:
        "0.75.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(
            [
                "results/{s.ID}/qc/samtools-stats/{s.sample}_processed.{s.genome_build}.txt",
                "results/{s.ID}/qc/fastqc/{s.sample}_processed.{s.genome_build}_fastqc.zip",
            ],
            s=list(samples.itertuples()),
        ),
    output:
        report(
            "results/{ID}/qc/multiqc.html",
            caption="workflow/report/multiqc.rst",
            category="Quality control",
        ),
    params:
        "--config config/multiqc_config.yaml"
    log:
        "results/logs/{ID}/multiqc/multiqc.log",
    wrapper:
        "0.75.0/bio/multiqc"
