rule fastqc:
    input:
        "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
    output:
        html="results/{ID}/qc/fastqc/{SAMPLE}.{GENOME}.html",
        zip="results/{ID}/qc/fastqc/{SAMPLE}.{GENOME}_fastqc.zip",
    log:
        "logs/{ID}/fastqc/fastqc_{SAMPLE}_all.{GENOME}.log",
    threads: 8
    resources:
        mem_mb = lambda wildcards,threads: 1024*threads
    wrapper:
        "v2.2.1/bio/fastqc"


rule samtools_stats:
    input:
        bam="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
    output:
        "results/{ID}/qc/samtools-stats/{SAMPLE}.{GENOME}.txt",
    log:
        "logs/{ID}/samtools_stats/samtools_stats_{SAMPLE}.{GENOME}.log",
    threads: 8
    wrapper:
        "v2.2.1/bio/samtools/stats"


rule mosdepth:
    input:
        bam="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
        bai="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
    output:
        "results/{ID}/qc/mosdepth/{SAMPLE}.{GENOME}.mosdepth.global.dist.txt",
        "results/{ID}/qc/mosdepth/{SAMPLE}.{GENOME}.mosdepth.region.dist.txt",
        "results/{ID}/qc/mosdepth/{SAMPLE}.{GENOME}.regions.bed.gz",
        summary="results/{ID}/qc/mosdepth/{SAMPLE}.{GENOME}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/{ID}/mosdepth/mosdepth_{SAMPLE}.{GENOME}.log",
    params:
        extra="--no-per-base",  # optional
        by="500",
    # additional decompression threads through `--threads`
    threads: 32  # This value - 1 will be sent to `--threads`
    wrapper:
        "v2.2.1/bio/mosdepth"


rule multiqc:
    input:
        expand(
            [
                "results/{s.ID}/qc/samtools-stats/{s.sample}.{s.genome_build}.txt",
                "results/{s.ID}/qc/fastqc/{s.sample}.{s.genome_build}_fastqc.zip",
                "results/{s.ID}/qc/mosdepth/{s.sample}.{s.genome_build}.mosdepth.global.dist.txt",
                "results/{s.ID}/qc/mosdepth/{s.sample}.{s.genome_build}.mosdepth.region.dist.txt",
            ],
            s=list(samples.itertuples()),
        ),
    output:
        report(
            "results/{ID}/qc/multiqc.html",
            caption="../report/QualityControl.rst",
            category="Quality control",
            labels={"Quality control": "MultiQC report"}
        ),
    log:
        "logs/{ID}/multiqc/multiqc.log",
    params:
        extra="--config config/multiqc_config.yaml",
        use_input_files_only=True,
    wrapper:
        "v2.6.0/bio/multiqc"
