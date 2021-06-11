
rule NGmerge:
    input:
        r1="results/fastq/{SAMPLE}_R1.fastq.gz",
        r2="results/fastq/{SAMPLE}_R2.fastq.gz",
    output:
        merged="results/NGmerge/{SAMPLE}_merged.fastq.gz",
        non_merged_1="results/NGmerge/{SAMPLE}_nonmerged_1.fastq.gz",
        non_merged_2="results/NGmerge/{SAMPLE}_nonmerged_2.fastq.gz",
    params:
        non_merged_prefix="results/NGmerge/{SAMPLE}_nonmerged",
    log:
        "results/logs/NGmerge/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        "workflow/scripts/NGmerge -n {threads} -z "
        "-1 {input.r1} -2 {input.r2} -o {output.merged} "
        "-f {params.non_merged_prefix} -l {log} "




rule map_reads:
    input:
        ref=lambda wc: config[wc.GENOME]["reference"],
        merged="results/NGmerge/{SAMPLE}_merged.fastq.gz",
        non_merged_1="results/NGmerge/{SAMPLE}_nonmerged_1.fastq.gz",
        non_merged_2="results/NGmerge/{SAMPLE}_nonmerged_2.fastq.gz",
    output:
        mapped_reads="results/mapped_reads/{SAMPLE}_all.{GENOME}.bam"
    params:
        RG=lambda wc: get_read_group(wc.SAMPLE),
    log:
        "results/logs/mapping/{SAMPLE}.{GENOME}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        "set +o pipefail;"
        "((bwa mem -t {threads} -R \"{params.RG}\" {input.ref} {input.merged}; "
        "bwa mem -t {threads} -R \"{params.RG}\" {input.ref} {input.non_merged_1} {input.non_merged_2} "
        "| grep -v \"^@\") | samtools view -Sbh -o {output.mapped_reads} - ) 2>{log}"

rule mark_duplicates:
    input:
        mapped_reads="results/mapped_reads/{SAMPLE}_all.{GENOME}.bam"
    output:
        processed_reads="results/mapped_reads/{SAMPLE}_processed.{GENOME}.bam"
    params:
        TMPdir=config["TMPdir"]
    log:"results/logs/markdup/{SAMPLE}.{GENOME}.log",
    conda:"../envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        "set +o pipefail;"
        "(samtools fixmate -u -m {input.mapped_reads} - | "
        "samtools sort -u -@ {threads} -T {params.TMPdir} - | "
        "samtools markdup -@ {threads} - {output.processed_reads}) 2>{log}"

rule index_bam:
    input:
        "results/mapped_reads/{SAMPLE}_processed.{GENOME}.bam"
    output:
        "results/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai"
    log:"results/logs/index_bam/{SAMPLE}.{GENOME}.log",
    conda: "../envs/cfDNA_prep.yaml"
    threads: 4
    shell:
        "samtools index -@ {threads} {input} {output} 2> {log}"