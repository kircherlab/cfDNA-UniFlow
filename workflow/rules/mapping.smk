rule map_reads_bwa2:
    input:
        unpack(get_mapping_input),
        ref=lambda wc: get_reference(wc),
        ref_index=lambda wc: multiext(get_reference(wc),".amb",".ann"),
    output:
        mapped_reads=temp("results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam")
    params:
        RG=lambda wc: get_read_group(wc.SAMPLE),
    log:
        "results/logs/{ID}/mapping/{SAMPLE}_all.{GENOME}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 32
    script:
        "../scripts/bwa-mem2_wrapper.py"

rule mark_duplicates:
    input:
        mapped_reads="results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam"
    output:
        processed_reads="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam"
    params:
        TMPDIR=config["TMPDIR"]
    log:"results/logs/{ID}/markdup/{SAMPLE}.{GENOME}.log",
    conda:"../envs/cfDNA_prep.yaml"
    threads: 64
    shell:
        #"set +o pipefail;"
        "(samtools fixmate -u -@ $(({threads}/3)) -m {input.mapped_reads} - | "
        "samtools sort -u -@ $(({threads}/3)) -T {params.TMPDIR} - | "
        "samtools markdup -@ $(({threads}/3)) - {output.processed_reads}) 2>{log}"

rule index_bam:
    input:
        "{path}.bam"
    output:
        "{path}.bam.bai"
    conda: "../envs/cfDNA_prep.yaml"
    threads: 32
    shell:
        "samtools index -@ {threads} {input} {output}"

