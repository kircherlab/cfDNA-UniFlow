rule map_reads_bwa2:
    input:
        unpack(get_mapping_input),
        ref=lambda wc: get_reference(wc),
        ref_index=lambda wc: multiext(get_reference(wc), ".amb", ".ann"),
    output:
        mapped_reads=temp("results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam"),
    params:
        RG=lambda wc: get_read_group(wc.SAMPLE),
    log:
        "logs/{ID}/map_reads_bwa2/map_reads_bwa2_{SAMPLE}_all.{GENOME}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 32
    script:
        "../scripts/bwa-mem2_wrapper.py"


rule mark_duplicates:
    input:
        unpack(get_mark_duplicates_input),
    output:
        processed_reads="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
    params:
        TMPDIR=config["TMPDIR"],
    log:
        "logs/{ID}/mark_duplicates/mark_duplicates_{SAMPLE}.{GENOME}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 64
    shell:
        "(samtools sort -n -@ $(({threads}/4)) -T {params.TMPDIR} {input.mapped_reads} | "
        "samtools fixmate -u -@ $(({threads}/4)) -m - - | "
        "samtools sort -u -@ $(({threads}/4)) -T {params.TMPDIR} - | "
        "samtools markdup -@ $(({threads}/4)) - {output.processed_reads}) 2>{log}"


rule index_bam:
    input:
        "{path}.bam",
    output:
        "{path}.bam.bai",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 32
    shell:
        "samtools index -@ {threads} {input} {output} 2>{log}"
