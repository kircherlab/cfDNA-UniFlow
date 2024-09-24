rule filter_merged:
    input:
        unfiltered="results/{ID}/NGmerge/merged/{SAMPLE}_merged.unfiltered.fastq.gz",
    output:
        filtered=temp("results/{ID}/NGmerge/merged/{SAMPLE}_merged.filtered.fastq.gz"),
    log:
        "logs/{ID}/filter_merged/filter_merged_{SAMPLE}.log",
    params:
        min_RL=config["length-filter"]["MINLEN"],
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        '(awk \'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}\' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}) 2> {log}'


rule filter_noadapter:
    input:
        unfiltered="results/{ID}/NGmerge/nonmerged/unfiltered.{SAMPLE}_noadapters_{read}.fastq.gz",
    output:
        filtered=temp(
            "results/{ID}/NGmerge/nonmerged/{SAMPLE}_noadapters_{read}.filtered.fastq.gz"
        ),
    log:
        "logs/{ID}/filter_noadapter/filter_noadapter_{SAMPLE}_noadapters_{read}.log",
    params:
        min_RL=config["length-filter"]["MINLEN"],
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        '(awk \'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}\' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}) 2> {log}'

rule filter_PE_singleton:
    input:
        unfiltered="results/{ID}/fastq/{SAMPLE}_PEsingleton.fastq.gz",
    output:
        filtered=temp("results/{ID}/fastq/{SAMPLE}_PEsingleton.filtered.fastq.gz"),
    log:
        "logs/{ID}/filter_PE_singleton/filter_PE_singleton_{SAMPLE}.log",
    params:
        min_RL=config["length-filter"]["MINLEN"],
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        '(awk \'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}\' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}) 2> {log}'

rule filter_single_read:
    input:
        unfiltered="results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.trimmed.fastq.gz",
    output:
        filtered=temp("results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.filtered.fastq.gz"),
    log:
        "logs/{ID}/filter_single_read/filter_single_read_{SAMPLE}.log",
    params:
        min_RL=config["length-filter"]["MINLEN"],
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        '(awk \'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}\' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}) 2> {log}'
