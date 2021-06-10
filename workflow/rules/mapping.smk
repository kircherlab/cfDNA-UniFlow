
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
        "scripts/NGmerge -n {threads} -z"
        "-1 {input.r1} -2 {input.r2} -o {output.merged}"
        "-f {params.non_merged_prefix} -l {log}"




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
        "bwa mem -t {threads} -R {params.RG} {input.ref} {input.merged};"
        "bwa mem -t {threads} -R {params.RG} {input.ref} {input.non_merged_1} {input.non_merged_2}"
        "| grep -v \"^@\" | samtools view -Sbh -o {output.mapped_reads} - ) 2>{log}"

