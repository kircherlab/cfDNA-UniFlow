
rule NGmerge:
    input:
        unpack(get_NGmerge_input),
        qual_table="resources/qual_profile.txt",
    output:
        merged=temp("results/{ID}/NGmerge/{SAMPLE}_merged.fastq.gz"),
        non_merged_1=temp("results/{ID}/NGmerge/{SAMPLE}_nonmerged_1.fastq.gz"),
        non_merged_2=temp("results/{ID}/NGmerge/{SAMPLE}_nonmerged_2.fastq.gz"),
    params:
        non_merged_prefix="results/{ID}/NGmerge/{SAMPLE}_nonmerged",
    log:
        "results/logs/{ID}/NGmerge/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        "set +o pipefail;"
        "workflow/scripts/NGmerge -w {input.qual_table} -u 41 -n {threads} -z "
        "-1 {input.r1} -2 {input.r2} -o {output.merged} "
        "-f {params.non_merged_prefix} -l {log} "




rule map_reads:
    input:
        ref=lambda wc: config[wc.GENOME]["reference"],
        merged="results/{ID}/NGmerge/{SAMPLE}_merged.fastq.gz",
        non_merged_1="results/{ID}/NGmerge/{SAMPLE}_nonmerged_1.fastq.gz",
        non_merged_2="results/{ID}/NGmerge/{SAMPLE}_nonmerged_2.fastq.gz",
    output:
        mapped_reads=temp("results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam")
    params:
        RG=lambda wc: get_read_group(wc.SAMPLE),
    log:
        "results/logs/{ID}/mapping/{SAMPLE}.{GENOME}.log",
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
        mapped_reads="results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam"
    output:
        processed_reads="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam"
    params:
        TMPdir=config["TMPdir"]
    log:"results/logs/{ID}/markdup/{SAMPLE}.{GENOME}.log",
    conda:"../envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        "set +o pipefail;"
        "(samtools fixmate -u -m {input.mapped_reads} - | "
        "samtools sort -u -@ {threads} -T {params.TMPdir} - | "
        "samtools markdup -@ {threads} - {output.processed_reads}) 2>{log}"

rule index_bam:
    input:
        "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam"
    output:
        config["output_root_path"] + "{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai"
    log:"results/logs/{ID}/index_bam/{SAMPLE}.{GENOME}.log",
    conda: "../envs/cfDNA_prep.yaml"
    threads: 4
    shell:
        "samtools index -@ {threads} {input} {output} 2> {log}"