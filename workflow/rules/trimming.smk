rule NGmerge:
    input:
        unpack(get_trimming_input),
        qual_table="resources/qual_profile.txt",
    output:
        merged=temp("results/{ID}/NGmerge/merged/{SAMPLE}_merged.unfiltered.fastq.gz"),
        non_merged_1=temp(
            "results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_1.fastq.gz"
        ),
        non_merged_2=temp(
            "results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_2.fastq.gz"
        ),
    params:
        non_merged_prefix="results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged",
        minlen=20,
    log:
        "logs/{ID}/NGmerge/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 64
    shell:
        "set +o pipefail;"
        "workflow/scripts/NGmerge -w {input.qual_table} -u 41 -d -e {params.minlen} -n {threads} -z "
        "-1 {input.r1} -2 {input.r2} -o {output.merged} "
        "-f {params.non_merged_prefix} -l {log} -v"


rule NGmerge_adapter:
    input:
        non_merged_1="results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_1.fastq.gz",
        non_merged_2="results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_2.fastq.gz",
        qual_table="resources/qual_profile.txt",
    output:
        R1_noadapters=temp(
            "results/{ID}/NGmerge/nonmerged/unfiltered.{SAMPLE}_noadapters_1.fastq.gz"
        ),
        R2_noadapters=temp(
            "results/{ID}/NGmerge/nonmerged/unfiltered.{SAMPLE}_noadapters_2.fastq.gz"
        ),
    params:
        adapt_minlen=1,
        output_prefix="results/{ID}/NGmerge/nonmerged/unfiltered.{SAMPLE}_noadapters",
    log:
        "logs/{ID}/NGmerge-adapter/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 64
    shell:
        "set +o pipefail;"
        "workflow/scripts/NGmerge -a -w {input.qual_table} -u 41 -d -e {params.adapt_minlen} -n {threads} -z "
        "-1 {input.non_merged_1} -2 {input.non_merged_2} -o {params.output_prefix}  "
        "-c {log} -v"


rule trimmomatic_pe:
    input:
        unpack(get_trimming_input),
    output:
        r1=temp("results/{ID}/trimmed/trimmomatic/{SAMPLE}.1.fastq.gz"),
        r2=temp("results/{ID}/trimmed/trimmomatic/{SAMPLE}.2.fastq.gz"),
        # reads where trimming entirely removed the mate
        r1_unpaired=temp(
            "results/{ID}/trimmed/trimmomatic/{SAMPLE}.1.unpaired.fastq.gz"
        ),
        r2_unpaired=temp(
            "results/{ID}/trimmed/trimmomatic/{SAMPLE}.2.unpaired.fastq.gz"
        ),
    log:
        "logs/{ID}/trimmomatic/{SAMPLE}.log",
    params:
        # list of trimmers (see manual)
        trimmer=get_trimmomatic_trimmers(),
        # optional parameters
        extra=get_trimmomatic_extra(),
        compression_level="-9",
    threads: 32
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v2.3.2/bio/trimmomatic/pe"


rule trimmomatic_se:
    input:
        unpack(get_SEtrimming_input)
    output:
        temp(
            "results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.trimmed.fastq.gz",
        ),
    log:
        "logs/{ID}/trimmomatic/{SAMPLE}.log",
    params:
        # list of trimmers (see manual)
        trimmer=get_trimmomatic_trimmers(),
        # optional parameters
        extra=get_trimmomatic_extra(),
        # optional compression levels from -0 to -9 and -11
        compression_level="-9",
    threads: 32
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v2.3.2/bio/trimmomatic/se"
