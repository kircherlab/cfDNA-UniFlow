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
        seed=config["SEED"],
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
        seed=config["SEED"],
        bl=lambda wildcards: "-bl {}".format(get_blacklist(wildcards)["blacklist"])
        if get_blacklist(wildcards)
        else "",
        normalized_interpolation="--normalized_interpolation"
        if config["GCbias_estimation"]["normalized_interpolation"]
        else "",
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
        {params.normalized_interpolation} \
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
            subcategory="GCbias-plot",
            labels={"Sample": "{SAMPLE}", "Type": "GCbias plot"},
        ),
    params:
        sample_name="{SAMPLE}",
        figsize="12 9",
    conda:
        "../envs/GC_bias.yaml"
    shell:
        """
        workflow/scripts/plot_GCbias_distribution.py \
        -o {output.GCbias_plot} \
        -s {params.sample_name} \
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
        -o {output.gc_weighted_bam} \
        -v
        """


def get_uncorrected_signals(wildcards):
    ID = wildcards.ID
    status_name = wildcards.status_name
    paths = expand(
        "-us results/{s.ID}/signals/signal-uncorrected/{target_region}.{s.sample}-uncorrected_{signal}.{s.genome_build}.csv.gz",
        s=list(
            samples.loc[
                (samples["ID"] == ID) & (samples["status"] == status_name)
            ].itertuples()
        ),
        allow_missing=True,
        target_region=wildcards.target_region,
        signal=wildcards.signal,
    )
    return paths


def get_corrected_signals(wildcards):
    ID = wildcards.ID
    status_name = wildcards.status_name
    paths = expand(
        "-cs results/{s.ID}/signals/signal-corrected/{target_region}.{s.sample}-corrected_{signal}.{s.genome_build}.csv.gz",
        s=list(
            samples.loc[
                (samples["ID"] == ID) & (samples["status"] == status_name)
            ].itertuples()
        ),
        target_region=wildcards.target_region,
        signal=wildcards.signal,
        allow_missing=True,
    )
    return paths


rule plot_GC_overlay:
    input:
        uncorrected_signals=lambda wc: expand(
            "results/{ID}/signals/signal-uncorrected/{target_region}.{s.sample}-uncorrected_{signal}.{s.genome_build}.csv.gz",
            s=list(
                samples.loc[
                    (samples["ID"] == wc.ID) & (samples["status"] == wc.status_name)
                ].itertuples()
            ),
            allow_missing=True,
        ),
        corrected_signals=lambda wc: expand(
            "results/{ID}/signals/signal-corrected/{target_region}.{s.sample}-corrected_{signal}.{s.genome_build}.csv.gz",
            s=list(
                samples.loc[
                    (samples["ID"] == wc.ID) & (samples["status"] == wc.status_name)
                ].itertuples()
            ),
            allow_missing=True,
        ),
    output:
        report(
            "results/{ID}/signals/GCcorrection-plots/{target_region}.{status_name}-GCcorrected_{signal}.{GENOME}.png",
            caption="../report/GCbias_overlay.rst",
            category="GCbias",
            subcategory="GCbias-region-overlay",
            labels={
                "Target region": "{target_region}",
                "Status": "{status_name}",
                "Type": "GCbias overlay",
            },
        ),
    params:
        uncorrected_samples=get_uncorrected_signals,
        corrected_samples=get_corrected_signals,
        name_regex="'.*\.(.*?)-[un]*corrected_.*'",  # matches everything between a '.' and '-[un]corrected_'
        signal="{signal}",
        target="{target_region}",
        overlay_mode=config["overlay_mode"],
        smoothing="--smoothing" if config["smoothing"] else "",
        smooth_window=config["smooth_window"],
        smooth_polyorder=config["smooth_polyorder"],
        rolling="--rolling" if config["rolling"] else "",
        rolling_window=config["rolling_window"],
        flank_norm="--flank_norm" if config["flank_norm"] else "",
        flank=config["flank"],
        display_window=config["display_window"],
        figsize=(12, 9),
        lower_limit="--lower_limit 0.8" if config["flank_norm"] else "",  # these options should be made configurable
        upper_limit="--upper_limit 1.2" if config["flank_norm"] else "",  # these options should be made configurable
    conda:
        "../envs/GC_bias.yaml"
    shell:
        """
        workflow/scripts/plot_overlay_GCcorrection.py \
        --output {output} \
        --regex {params.name_regex} \
        --signal {params.signal} \
        --target {params.target} \
        --overlay_mode {params.overlay_mode} \
        {params.smoothing} \
        --smooth_window {params.smooth_window} \
        --smooth_polyorder {params.smooth_polyorder} \
        {params.rolling} \
        --rolling_window {params.rolling_window} \
        {params.flank_norm} \
        --flank {params.flank} \
        --display_window {params.display_window} \
        --figsize {params.figsize} \
        {params.lower_limit} \
        {params.upper_limit} \
        {params.uncorrected_samples} \
        {params.corrected_samples}
        """
