
#! wrapper rules können auch anders benannte inputs und outputs haben -> output based on configs!

#rule get_genome:
#    output:
#        "resources/genome.fasta",
#    log:
#        "logs/get-genome.log",
#    params:
#        species=config["ref"]["species"],
#        datatype="dna",
#        build=config["ref"]["build"],
#        release=config["ref"]["release"],
#    cache: True
#    wrapper:
#        "0.75.0/bio/reference/ensembl-sequence"

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