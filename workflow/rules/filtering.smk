rule filter_merged:
    input:
        unfiltered="results/{ID}/NGmerge/merged/{SAMPLE}_merged.unfiltered.fastq.gz",
    output:
        filtered=temp("results/{ID}/NGmerge/merged/{SAMPLE}_merged.filtered.fastq.gz"),
    params:
        min_RL=config["length-filter"]["MINLEN"],
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        'awk \'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}\' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}'


rule filter_noadapter:
    input:
        unfiltered="results/{ID}/NGmerge/nonmerged/unfiltered.{SAMPLE}_noadapters_{read}.fastq.gz",
    output:
        filtered=temp(
            "results/{ID}/NGmerge/nonmerged/{SAMPLE}_noadapters_{read}.filtered.fastq.gz"
        ),
    params:
        min_RL=config["length-filter"]["MINLEN"],
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        'awk \'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}\' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}'

rule filter_PE_singleton:
    input:
        unfiltered="results/{ID}/fastq/{SAMPLE}_PEsingleton.fastq.gz",
    output:
        filtered=temp("results/{ID}/fastq/{SAMPLE}_PEsingleton.filtered.fastq.gz"),
    params:
        min_RL=config["length-filter"]["MINLEN"],
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        'awk \'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}\' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}'

rule filter_single_read:
    input:
        unfiltered="results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.trimmed.fastq.gz",
    output:
        filtered=temp("results/{ID}/trimmed/trimmomatic/{SAMPLE}_single_read.filtered.fastq.gz"),
    params:
        min_RL=config["length-filter"]["MINLEN"],
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        'awk \'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}\' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}'
