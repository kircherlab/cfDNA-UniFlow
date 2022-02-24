
#! wrapper rules können auch anders benannte inputs und outputs haben -> output based on configs!


rule get_reference:
    output:
        "resources/reference/{GENOME}.fa",
    log:
        "logs/get-{GENOME}-reference.log",
    params:
        url= lambda wc:get_ref_url(wc)
    shell:
        "(curl -L {params.url} | gzip -d > {output}) 2> {log}"
        


rule bwa_index:
    input:
        ref="resources/reference/{GENOME}.fa",
    output:
        ref=multiext( "resources/reference/{GENOME}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        algorithm="bwtsw",
    log:
       "logs/bwa_index-{GENOME}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    threads: 1
    shell:
        "bwa index -a {params.algorithm} {input.ref}"
