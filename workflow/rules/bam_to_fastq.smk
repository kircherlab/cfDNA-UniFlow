# Convert bam file to R1 and R2 fastQ files.
# Singletons (reads not having a partner read) will be saved to a separate file.

rule bam_to_fastq:
    input:
        bam=lambda wildcards: samples["bam"][wildcards.SAMPLE],
        bai=lambda wildcards: samples["bam"][wildcards.SAMPLE]+".bai",
    output:
        r1=temp("results/{ID}/fastq/{SAMPLE}_R1.fastq.gz"),
        r2=temp("results/{ID}/fastq/{SAMPLE}_R2.fastq.gz"),
        s1=temp("results/{ID}/fastq/{SAMPLE}_single_read.fastq.gz"),
    log:
        "results/logs/{ID}/bam_to_fastq/{SAMPLE}.log",
    params:
        TMPDIR=config["TMPDIR"],
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 32
    shell:
        """set +o pipefail;
        ( samtools collate -T {params.TMPDIR} -@ $(({threads}/2)) -O {input.bam} | \
        samtools fastq -@ $(({threads}/2)) -t \
        -1 {output.r1} \
        -2 {output.r2} \
        -s {output.s1} \
        )2> {log}"""
