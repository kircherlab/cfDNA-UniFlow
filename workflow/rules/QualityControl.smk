
rule fastqc:
    input:
        "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
    output:
        html="results/{ID}/qc/fastqc/{SAMPLE}_all.{GENOME}.html",
        zip="results/{ID}/qc/fastqc/{SAMPLE}_all.{GENOME}_fastqc.zip",
    log:
        "results/logs/{ID}/fastqc/{SAMPLE}_all.{GENOME}.log",
    wrapper:
        "0.75.0/bio/fastqc"


rule samtools_stats:
    input:
        "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
    output:
        "results/{ID}/qc/samtools-stats/{SAMPLE}_all.{GENOME}.txt",
    log:
        "results/logs/{ID}/fastqc/{SAMPLE}.{GENOME}.log",
    wrapper:
        "0.75.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(
            [
                "results/{s.ID}/qc/samtools-stats/{s.sample}_all.{s.genome_build}.txt",
                "results/{s.ID}/qc/fastqc/{s.sample}_all.{s.genome_build}_fastqc.zip",
            ],
            s=list(samples.itertuples()),
        ),
    output:
        report(
            "results/{ID}/qc/multiqc.html",
            caption="workflow/report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "results/logs/{ID}/multiqc/multiqc.log",
    wrapper:
        "0.75.0/bio/multiqc"
