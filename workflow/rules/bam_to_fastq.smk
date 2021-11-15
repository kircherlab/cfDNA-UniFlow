ruleorder: flagstat_pre_conversion>bam_to_fastq

rule flagstat_pre_conversion:
    input:
        bam=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        "results/{ID}/flagstat/{SAMPLE}_flagstat.txt.gz"
    log:
        "results/logs/{ID}/flagstat/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 16
    shell:
        "set +o pipefail;"
        "samtools flagstat -@ {threads} {input.bam} 2> {log} | gzip -c > {output}"


rule bam_to_fastq:
    input:
        bam=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        r1=temp("results/{ID}/fastq/{SAMPLE}_R1.fastq.gz"),
        r2=temp("results/{ID}/fastq/{SAMPLE}_R2.fastq.gz"),
        s1=temp("results/{ID}/fastq/{SAMPLE}_single_read.fastq.gz"),
    log:
        "results/logs/{ID}/bam_to_fastq/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 16
    shell:
        """set +o pipefail;
        (
                    samtools sort -n -@ {threads} {input.bam} | \
        samtools fastq -@ {threads} -t \
        -1 {output.r1} \
        -2 {output.r2} \
        -0 {output.s1} \
        )2> {log}"""