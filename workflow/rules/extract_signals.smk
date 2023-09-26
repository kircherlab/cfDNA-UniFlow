

rule extract_counts:
    input:
        #target="results/regions/{GENOME}/target_region/{target_region}.blacklist-excluded.bed.gz",
        target=lambda wildcards: regions["path"][wildcards.target_region],
        BAMFILE="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
        BAIFILE="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
        twobit_genome=lambda wildcards: config[wildcards.GENOME]["2bit_ref"],
    output:
        WPS="results/{ID}/signals/signal-uncorrected/{target_region}.{SAMPLE}-uncorrected_WPS.{GENOME}.csv.gz",
        COV="results/{ID}/signals/signal-uncorrected/{target_region}.{SAMPLE}-uncorrected_COV.{GENOME}.csv.gz",
        MEAN_WEIGHT="results/{ID}/signals/signal-uncorrected/{target_region}.{SAMPLE}-uncorrected_MEAN_WEIGHT.{GENOME}.csv.gz",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        lengthSR=76,
        protection=120,
        seed="--seed 42",
        out_pre="results/{ID}/signals/signal-uncorrected/{target_region}.{SAMPLE}-uncorrected_%s.{GENOME}.csv.gz",
    threads: 30
    conda:
        #"../envs/cfDNA_rework.yml"
        "../envs/GC_bias.yaml"
    shell:
        """
        workflow/scripts/extract_cfDNA_signals.py \
        -b {input.BAMFILE} \
        -g {input.twobit_genome} \
        -r {input.target} \
        -o {params.out_pre} \
        -min {params.minRL} \
        -max {params.maxRL} \
        -w {params.protection} \
        -l {params.lengthSR} \
        -p {threads} \
        {params.seed}
        """


rule extract_GCcorrected_counts:
    input:
        #target="results/regions/{GENOME}/target_region/{target_region}.blacklist-excluded.bed.gz",
        target=lambda wildcards: regions["path"][wildcards.target_region],
        BAMFILE="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
        BAIFILE="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
        twobit_genome=lambda wildcards: config[wildcards.GENOME]["2bit_ref"],
        #GCfreqfile="results/{ID}/GCBias/bias_table/{SAMPLE}-{computation}_GCbias_interpolated_{analysis}.{GENOME}.tsv.gz",
    output:
        WPS="results/{ID}/signals/signal-corrected/{target_region}.{SAMPLE}-corrected_WPS.{GENOME}.csv.gz",
        COV="results/{ID}/signals/signal-corrected/{target_region}.{SAMPLE}-corrected_COV.{GENOME}.csv.gz",
        MEAN_WEIGHT="results/{ID}/signals/signal-corrected/{target_region}.{SAMPLE}-corrected_MEAN_WEIGHT.{GENOME}.csv.gz",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        lengthSR=config["lengthSR"],
        protection=config["bpProtection"],
        seed="--seed 42",
        out_pre="results/{ID}/signals/signal-corrected/{target_region}.{SAMPLE}-corrected_%s.{GENOME}.csv.gz",
    threads: 30
    conda:
        "../envs/GC_bias.yaml"
    shell:
        """
        workflow/scripts/extract_cfDNA_signals.py \
        -b {input.BAMFILE} \
        -g {input.twobit_genome} \
        -r {input.target} \
        -o {params.out_pre} \
        -gc \
        -min {params.minRL} \
        -max {params.maxRL} \
        -w {params.protection} \
        -l {params.lengthSR} \
        -p {threads} \
        {params.seed}
        """


def get_control_files(wildcards):
    control_name = config["control_name"]
    path_template = "results/{s.ID}/signals/signal-{correction}/{target_region}.{s.sample}-{correction}_{signal}.{s.genome_build}.csv.gz"
    paths =  expand(
        path_template,
        s=list(samples.loc[samples["status"] == control_name].itertuples()),
        allow_missing=True,
    )
    return paths


rule aggregate_controls:
    input:
        control_files=get_control_files,
        #control_files=["results/cfDNA_UniFlow/signals/signal-corrected/LYL1.EGAF00002727162-corrected_COV.hg38.csv.gz",
        #        "results/cfDNA_UniFlow/signals/signal-corrected/LYL1.EGAF00002727163-corrected_COV.hg38.csv.gz"]
    output:
        "results/{ID}/signals/controls/{target_region}.{control_name}-{correction}_{signal}.{GENOME}.tsv.gz",
    params:
        control_name=config["control_name"],
        overlay_mode=config["overlay_mode"],
        smoothing= "--smoothing" if config["smoothing"] else "",
        smooth_window=config["smooth_window"],
        smooth_polyorder=config["smooth_polyorder"],
        rolling="--rolling" if config["rolling"] else "",
        rolling_window=config["rolling_window"],
        flank_norm="--flank_norm" if config["flank_norm"] else "",
        flank=config["flank"],
    conda:
        "../envs/GC_bias.yaml"
    shell:
        """
        workflow/scripts/aggregate_controls.py \
        --output {output}\
        --control_name {params.control_name} \
        --overlay_mode {params.overlay_mode} \
        {params.smoothing} \
        --smooth_window {params.smooth_window} \
        --smooth_polyorder {params.smooth_polyorder} \
        {params.rolling} \
        --rolling_window {params.rolling_window} \
        {params.flank_norm} \
        --flank {params.flank} \
        {input.control_files}
        """

def get_case_files(wildcards):
    case_name = wildcards.case_name
    path_template = "results/{s.ID}/signals/signal-{correction}/{target_region}.{s.sample}-{correction}_{signal}.{s.genome_build}.csv.gz"
    paths =  expand(
        path_template,
        s=list(samples.loc[samples["status"] == case_name].itertuples()),
        allow_missing=True,
    )
    return paths

rule plot_case_control:
    input:
        control_file= "results/{ID}/signals/controls/{target_region}.{control_name}-{correction}_{signal}.{GENOME}.tsv.gz",
        case_files = get_case_files,
    output:
            CaseControl_plot=report(
            "results/{ID}/signals/case-control-plots/{target_region}.{case_name}-vs-{control_name}-{correction}_{signal}.{GENOME}.png",
            caption="../report/Case-control_plot.rst",
            category="Case-control",
            labels={"Target region": "{target_region}", "Case": "{case_name}", "Control":"{control_name}", "Type": "Case-control overlay"},
        ),
    params:
        signal = "{signal}",
        target = "{target_region}",
        case_name="{case_name}",
        overlay_mode=config["overlay_mode"],
        smoothing= "--smoothing" if config["smoothing"] else "",
        smooth_window=config["smooth_window"],
        smooth_polyorder=config["smooth_polyorder"],
        rolling="--rolling" if config["rolling"] else "",
        rolling_window=config["rolling_window"],
        flank_norm="--flank_norm" if config["flank_norm"] else "",
        flank=config["flank"],
        display_window = config["display_window"],
        aggregate_controls = "--aggregate_controls" if config["aggregate_controls"] else "",
        GC_corrected = "--GC_corrected" if config["utility"]["GCbias-correction"] else "",
        figsize = (12, 9),
        lower_limit="--lower_limit 0.8" if config["flank_norm"] else "", # these options should be made configurable
        upper_limit="--upper_limit 1.2" if config["flank_norm"] else "", # these options should be made configurable
    conda:
        "../envs/GC_bias.yaml"
    shell:
        """
        workflow/scripts/plot_overlay_case-control.py \
        --control_samples {input.control_file} \
        --output {output.CaseControl_plot} \
        --signal {params.signal} \
        --target {params.target} \
        --case_name {params.case_name} \
        --overlay_mode {params.overlay_mode} \
        {params.smoothing} \
        --smooth_window {params.smooth_window} \
        --smooth_polyorder {params.smooth_polyorder} \
        {params.rolling} \
        --rolling_window {params.rolling_window} \
        {params.flank_norm} \
        --flank {params.flank} \
        --display_window {params.display_window} \
        {params.aggregate_controls} \
        {params.GC_corrected} \
        --figsize {params.figsize} \
        {params.lower_limit} \
        {params.upper_limit} \
        {input.case_files}
        """

