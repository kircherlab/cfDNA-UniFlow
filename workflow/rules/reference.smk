
rule get_reference:
    output:
        "resources/reference/{GENOME}.fa",
    log:
        "logs/get-{GENOME}-reference.log",
    params:
        url= lambda wc:get_ref_url(wc)
    shell:
        "(curl -L {params.url} | gzip -d > {output}) 2> {log}"
        


rule bwa_mem2_index:
    input:
        ref="resources/reference/{GENOME}.fa",
    output:
        ref=multiext( "resources/reference/{GENOME}.fa", ".amb", ".ann", ".pac", ),
        #ref=multiext( "resources/reference/{GENOME}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
       "logs/bwa_index-{GENOME}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 16
    shell:
        "bwa-mem2 index {input.ref}"


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
        URLs = [
        "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/NexteraPE-PE.fa",
        "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq2-PE.fa",
        "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq2-SE.fa",
        "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE-2.fa",
        "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa",
        "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-SE.fa",
        ]
    log:
        "logs/get_trimmomatic_adapter.log"
    shell:
        "wget -P {params.prefix} {params.URLs} -o {log}"
