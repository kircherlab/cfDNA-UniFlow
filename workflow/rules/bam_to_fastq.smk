# Convert bam file to R1 and R2 fastQ files.
# Singletons (reads not having a partner read) will be saved to a separate file.


rule bam_to_fastq:
    input:
        bam=lambda wildcards: samples["bam"][wildcards.SAMPLE],
        bai=lambda wildcards: samples["bam"][wildcards.SAMPLE] + ".bai",
    output:
        r1=temp("results/{ID}/fastq/{SAMPLE}_R1.fastq.gz"), # Read 1 in Paired End Sequencing
        r2=temp("results/{ID}/fastq/{SAMPLE}_R2.fastq.gz"), # Read 2 in Paired End Sequencing
        s1=temp("results/{ID}/fastq/{SAMPLE}_PEsingleton.fastq.gz"), # Singletons (paired, but only one read present) in Paired End Sequencing
        SR=temp("results/{ID}/fastq/{SAMPLE}_single_read.fastq.gz"), # File for Single End Sequencing (flags for both R1 and R2 set, or both unset)
    log:
        "logs/{ID}/bam_to_fastq/bam_to_fastq_{SAMPLE}.log",
    params:
        TMPDIR=config["TMPDIR"]+"{SAMPLE}",
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
        -0 {output.SR} \
        )2> {log}"""
