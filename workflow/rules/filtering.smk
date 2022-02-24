rule filter_merged:
    input:
        unfiltered="results/{ID}/NGmerge/merged/{SAMPLE}_merged.unfiltered.fastq.gz"
    output:
        filtered=temp("results/{ID}/NGmerge/merged/{SAMPLE}_merged.filtered.fastq.gz")
    params:
        min_RL = config["NGmerge"]["MINLEN"]
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        "awk 'BEGIN {{FS = \"\\t\" ; OFS = \"\\n\"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}"


rule filter_interleaved:
    input:
        unfiltered="results/{ID}/NGmerge/nonmerged/{SAMPLE}_interleaved_noadapters.unfiltered.fastq.gz"
    output:
        filtered=temp("results/{ID}/NGmerge/nonmerged/{SAMPLE}_interleaved_noadapters.filtered.fastq.gz"),
    params:
        min_RL = config["NGmerge"]["MINLEN"]
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        "awk 'BEGIN {{FS = \"\\t\" ; OFS = \"\\n\"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}"

rule filter_single_read:
    input:
        unfiltered="results/{ID}/fastq/{SAMPLE}_single_read.fastq.gz"
    output:
        filtered=temp("results/{ID}/fastq/{SAMPLE}_single_read.filtered.fastq.gz")
    params:
        min_RL = config["NGmerge"]["MINLEN"]
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        "awk 'BEGIN {{FS = \"\\t\" ; OFS = \"\\n\"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}"