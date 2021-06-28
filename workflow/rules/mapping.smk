
rule NGmerge:
    input:
        unpack(get_NGmerge_input),
        qual_table="resources/qual_profile.txt",
    output:
        merged="results/{ID}/NGmerge/merged/{SAMPLE}_merged.unfiltered.fastq.gz",#temp("results/{ID}/NGmerge/merged/{SAMPLE}_merged.fastq.gz"),
        non_merged_1="results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_1.fastq.gz",#temp("results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_1.fastq.gz"),
        non_merged_2="results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_2.fastq.gz",#temp("results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_2.fastq.gz"),
    params:
        non_merged_prefix="results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged",
        minlen=20
    log:
        "results/logs/{ID}/NGmerge/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 8
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
        interleaved_output="results/{ID}/NGmerge/nonmerged/{SAMPLE}_interleaved_noadapters.unfiltered.fastq.gz"
        #noadapters_1="results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_noadapters_1.fastq.unfiltered.gz",
        #noadapters_2="results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_noadapters_2.fastq.unfiltered.gz"#temp("results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_noadapters.fastq.gz"),
    params:
        #noadapter_prefix="results/{ID}/NGmerge/nonmerged/{SAMPLE}_nonmerged_noadapters",
        adapt_minlen=1
    log:
        "results/logs/{ID}/NGmerge-adapter/{SAMPLE}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        "set +o pipefail;"
        "workflow/scripts/NGmerge -a -i -w {input.qual_table} -u 41 -d -e {params.adapt_minlen} -n {threads} -z "
        "-1 {input.non_merged_1} -2 {input.non_merged_2} -o {output.interleaved_output}  "
        "-c {log} -v"

rule filter_interleaved:
    input:
        unfiltered="results/{ID}/NGmerge/nonmerged/{SAMPLE}_interleaved_noadapters.unfiltered.fastq.gz"
    output:
        filtered="results/{ID}/NGmerge/nonmerged/{SAMPLE}_interleaved_noadapters.filtered.fastq.gz"
    params:
        min_RL = 30
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        "awk 'BEGIN {{FS = \"\\t\" ; OFS = \"\\n\"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}"

rule filter_merged:
    input:
        unfiltered="results/{ID}/NGmerge/merged/{SAMPLE}_merged.unfiltered.fastq.gz"
    output:
        filtered="results/{ID}/NGmerge/merged/{SAMPLE}_merged.filtered.fastq.gz"
    params:
        min_RL = 30
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        "awk 'BEGIN {{FS = \"\\t\" ; OFS = \"\\n\"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.min_RL}) {{print header, seq, qheader, qseq}}}}' <( zcat {input.unfiltered} ) |  gzip -c > {output.filtered}"



rule map_reads:
    input:
        ref=lambda wc: config[wc.GENOME]["reference"],
        merged="results/{ID}/NGmerge/{SAMPLE}_merged.fastq.gz",
        non_merged_1="results/{ID}/NGmerge/{SAMPLE}_nonmerged_1.fastq.gz",
        non_merged_2="results/{ID}/NGmerge/{SAMPLE}_nonmerged_2.fastq.gz",
    output:
        mapped_reads=temp("results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam")
    params:
        RG=lambda wc: get_read_group(wc.SAMPLE),
    log:
        "results/logs/{ID}/mapping/{SAMPLE}.{GENOME}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        "set +o pipefail;"
        "((bwa mem -t {threads} -R \"{params.RG}\" {input.ref} {input.merged}; "
        "bwa mem -t {threads} -R \"{params.RG}\" {input.ref} {input.non_merged_1} {input.non_merged_2} "
        "| grep -v \"^@\") | samtools view -Sbh -o {output.mapped_reads} - ) 2>{log}"

rule mark_duplicates:
    input:
        mapped_reads="results/{ID}/mapped_reads/{SAMPLE}_all.{GENOME}.bam"
    output:
        processed_reads="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam"
    params:
        TMPdir=config["TMPdir"]
    log:"results/logs/{ID}/markdup/{SAMPLE}.{GENOME}.log",
    conda:"../envs/cfDNA_prep.yaml"
    threads: 8
    shell:
        "set +o pipefail;"
        "(samtools fixmate -u -m {input.mapped_reads} - | "
        "samtools sort -u -@ {threads} -T {params.TMPdir} - | "
        "samtools markdup -@ {threads} - {output.processed_reads}) 2>{log}"

rule index_bam:
    input:
        "results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam"
    output:
        config["output_root_path"] + "{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai"
    log:"results/logs/{ID}/index_bam/{SAMPLE}.{GENOME}.log",
    conda: "../envs/cfDNA_prep.yaml"
    threads: 4
    shell:
        "samtools index -@ {threads} {input} {output} 2> {log}"