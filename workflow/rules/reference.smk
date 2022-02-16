
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
        lambda wc: config[wc.GENOME]["reference"],
    output:
        ref=multiext(config[wc.GENOME]["reference"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
        cache: True
    wrapper:
        "0.74.0/bio/bwa/index"