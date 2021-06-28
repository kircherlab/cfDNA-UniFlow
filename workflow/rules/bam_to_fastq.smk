
rule bam_to_fastq:
    input:
        bam=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        r1=temp("results/{ID}/fastq/{SAMPLE}_R1.fastq.gz"),
        r2=temp("results/{ID}/fastq/{SAMPLE}_R2.fastq.gz"),
        s1=temp("results/{ID}/fastq/{SAMPLE}_singleton.fastq.gz"),
    log:
        "results/logs/{ID}/bam_to_fastq/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        """(
            samtools sort -n -@ {threads} {input.bam} | samtools fastq -@ {threads} -t -1 {output.r1} -2 {output.r2} -s {output.s1}
        )2> {log}
        """
