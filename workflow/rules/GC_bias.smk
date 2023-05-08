def get_blacklist(wildcards):
    if wildcards.blacklist == "repeatmasker":
        file = config[wildcards.GENOME]["repeatmasker"]
        return {"blacklist": file}
    else:
        return {}


def get_ref_bam_whole_genome(wildcards):
    """Gets a reference bam file for a given genome build.

    Args:
        wildcards (_type_): _description_

    Returns:
        _type_: _description_
    """

    sample = samples["sample"].loc[samples.genome_build == wildcards.GENOME].iloc[0]
    bamfile = "results/{{ID}}/mapped_reads/{SAMPLE}_all.{{GENOME}}.bam".format(
        SAMPLE=sample
    )
    baifile = bamfile + ".bai"
    return {"BAMFILE": bamfile, "BAIFILE": baifile}


rule computeGCbias_background:
    input:
        unpack(get_blacklist),
        unpack(get_ref_bam_whole_genome),
        twobit_genome=lambda wildcards: config[wildcards.GENOME]["2bit_ref"],
    output:
        background="results/{ID}/GCBias_background/precomputed_background_{blacklist}.{GENOME}.tsv.gz",
    params:
        seed=42,
        bl=lambda wildcards: "-bl {}".format(get_blacklist(wildcards)["blacklist"])
        if get_blacklist(wildcards)
        else "",
    log:
        "logs/{ID}/computeGCbias_background/precomputed_background_{blacklist}.{GENOME}.log",
    conda:
        "../envs/GC_bias.yaml"
    threads: 48
    shell:
        """
        computeGCBias_background -b {input.BAMFILE} \
        -g {input.twobit_genome} \
        -i \
        {params.bl} \
        -p {threads} \
        -mp=MPIRE \
        --output {output.background} \
        --seed {params.seed} \
        --standard_chroms \
        -v
        """


rule computeGCbias:
    input:
        unpack(get_blacklist),
        BAMFILE="results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam",
        BAIFILE="results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam" + ".bai",
        background="results/{ID}/GCBias_background/precomputed_background_{blacklist}.{GENOME}.tsv.gz",
        twobit_genome=lambda wildcards: config[wildcards.GENOME]["2bit_ref"],
    output:
        GCfreqfile="results/{ID}/GCBias/bias_table/{SAMPLE}-GCbias_{blacklist}.{GENOME}.tsv.gz",
    params:
        seed=42,
        bl=lambda wildcards: "-bl {}".format(get_blacklist(wildcards)["blacklist"])
        if get_blacklist(wildcards)
        else "",
        #read_length = 150
    log:
        "logs/{ID}/computeGCbias/{SAMPLE}-GCbias_{blacklist}.{GENOME}.log",
    conda:
        "../envs/GC_bias.yaml"
    threads: 48
    benchmark:
        "results/benchmarks/{ID}/benchmark_whole_genome/computeGCbias_interpolated_CSAPS_precomputed_background_{blacklist}-{GENOME}_{SAMPLE}-resources.tsv"
    shell:
        """
        computeGCBias_readlen -b {input.BAMFILE} \
        -g {input.twobit_genome} \
        -i \
        {params.bl} \
        --precomputed_background {input.background} \
        -p {threads} \
        -mp=MPIRE \
        --GCbiasFrequenciesFile {output.GCfreqfile} \
        --seed {params.seed} \
        --standard_chroms \
        -v &> {log}
        """

rule plot_GCbias:
    input:
        GCfreqfile="results/{ID}/GCBias/bias_table/{SAMPLE}-GCbias_{blacklist}.{GENOME}.tsv.gz",
    output:
        GCbias_plot=report(
            "results/{ID}/GCBias/plots/{SAMPLE}-GCbias-plot_{blacklist}.{GENOME}.png",
            caption="../report/GCbias.rst",
            category="GCbias",
            labels={
              "Sample": "{SAMPLE}",
              "Type":"GCbias plot"
          }
        ),
    params:
        sample_name="{SAMPLE}",
        quantile_threshold="0.999",
        figsize="15 12",
    conda:
        "../envs/GC_bias.yaml"
    shell:
        """
        workflow/scripts/plot_GCbias_distribution.py \
        -o {output.GCbias_plot} \
        -s {params.sample_name} \
        -q {params.quantile_threshold} \
        --figsize {params.figsize} \
        {input.GCfreqfile}
        """


rule correctGCbias:
    input:
        BAMFILE="results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam",
        BAIFILE="results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam" + ".bai",
        twobit_genome=lambda wildcards: config[wildcards.GENOME]["2bit_ref"],
        GCfreqfile="results/{ID}/GCBias/bias_table/{SAMPLE}-GCbias_{blacklist}.{GENOME}.tsv.gz",
    output:
        gc_weighted_bam=temp(
            "results/{ID}/corrected_reads/{SAMPLE}_GCcorrected_{blacklist}.{GENOME}.bam"
        ),
        gc_weighted_bai=temp(
            "results/{ID}/corrected_reads/{SAMPLE}_GCcorrected_{blacklist}.{GENOME}.bam.bai"
        ),
    params:
        GC_weights="-w",
        effectiveGenomeSize=lambda wildcards: config[wildcards.GENOME][
            "effectiveGenomeSize"
        ],
    log:
        "logs/{ID}/correctGCbias/{SAMPLE}-GCbias_{blacklist}.{GENOME}.log",
    conda:
        "../envs/GC_bias.yaml"
    threads: 48
    shell:
        """
        correctGCBias_readlen -b {input.BAMFILE} \
        --effectiveGenomeSize {params.effectiveGenomeSize} \
        -g {input.twobit_genome} \
        --GCbiasFrequenciesFile {input.GCfreqfile} \
        -p {threads} \
        {params.GC_weights} \
        -o {output.gc_weighted_bam}
        """
