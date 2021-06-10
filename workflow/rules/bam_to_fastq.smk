
rule bam_to_fastq:
    input:
        bam=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        r1="results/fastq/{SAMPLE}_R1.fastq.gz",
        r2="results/fastq/{SAMPLE}_R2.fastq.gz",
    log:
        "results/logs/bam_to_fastq/{SAMPLE}.log",
    conda:
        "workflow/envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        """(
            samtools sort -n -@ {threads} {input.bam} | samtools fastq -@ {threads} -t -1 {output.r1} -2 {output.r2}
        )2> {log}
        """
