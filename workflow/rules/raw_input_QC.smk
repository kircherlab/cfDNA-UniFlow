ruleorder: flagstat_pre_conversion>bam_to_fastq


## Flagstat for bam files
rule flagstat_pre_conversion:
    input:
        bam=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        "results/{ID}/flagstat/{SAMPLE}_flagstat.txt.gz"
    log:
        "results/logs/{ID}/flagstat/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 64
    shell:
        "set +o pipefail;"
        "samtools flagstat -@ {threads} {input.bam} 2> {log} | gzip -c > {output}"

## samtools stats for bam files


## FastQC: BAM, SAM or FastQ files