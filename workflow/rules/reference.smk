
rule get_fasta_reference:
    output:
        "resources/reference/{GENOME}.fa",
    log:
        "logs/get_fasta_reference_{GENOME}.log",
    params:
        url=lambda wc: get_fasta_ref_url(wc),
        tmp_name="resources/reference/{GENOME}.fa.gz",
    shell:
        "(rsync -avzP {params.url} {params.tmp_name}; gzip -f -d {params.tmp_name} > {output} )2> {log}"


rule get_twobit_reference:
    output:
        "resources/reference/{GENOME}.2bit",
    log:
        "logs/get_twobit_reference_{GENOME}.log",
    params:
        url=lambda wc: get_twobit_ref_url(wc),
    shell:
        "rsync -avzP {params.url} {output} 2> {log}"


rule bwa_mem2_index:
    input:
        ref="resources/reference/{GENOME}.fa",
    output:
        ref=multiext(
            "resources/reference/{GENOME}.fa",
            ".amb",
            ".ann",
            ".pac",
        ),
    log:
        "logs/bwa_mem2_index-{GENOME}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 16
    shell:
        "bwa-mem2 index {input.ref} 2> {log}"


rule get_trimmomatic_adapters:
    output:
        "resources/adapter/NexteraPE-PE.fa",
        "resources/adapter/TruSeq2-PE.fa",
        "resources/adapter/TruSeq2-SE.fa",
        "resources/adapter/TruSeq3-PE-2.fa",
        "resources/adapter/TruSeq3-PE.fa",
        "resources/adapter/TruSeq3-SE.fa",
    params:
        prefix="resources/adapter/",
        URLs=[
            "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/NexteraPE-PE.fa",
            "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq2-PE.fa",
            "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq2-SE.fa",
            "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE-2.fa",
            "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa",
            "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-SE.fa",
        ],
    log:
        "logs/get_trimmomatic_adapter.log",
    shell:
        "wget -P {params.prefix} {params.URLs} -o {log}"


rule get_ichorCNA_files:
    output:
        zip="ichorCNA_master.zip",
        dir=directory("resources/ichorCNA"),
    params:
        prefix="resources/",
        URL="https://github.com/broadinstitute/ichorCNA/archive/refs/heads/master.zip",  #"https://github.com/GavinHaLab/ichorCNA/archive/refs/heads/master.zip"
    log:
        "logs/get_ichorCNA_files.log",
    shell:
        '(wget -P {params.prefix} -c {params.URL} -O {output.zip}; unzip -j {output.zip} "ichorCNA-master/inst/*" -d {output.dir}) 2> {log}'
