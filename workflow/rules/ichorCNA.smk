
rule get_chroms:
    input:
        bam="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
        bai="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
    output:
        chroms="results/{ID}/icorCNA/chroms/{SAMPLE}_processed.{GENOME}.chromosomes.txt",
    log:
        "results/logs/{ID}/get_chroms/{SAMPLE}_processed.{GENOME}.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        """samtools idxstats {input.bam} | \
        cut -f 1 | grep -w 'chr[1-9]\|chr[1-2][0-9]\|chr[X,Y]\|^[1-2][0-9]\|^[0-9]\|^[X,Y]' \
        | tr '\n' ','| sed 's/,*\r*$//' \
        1> {output.chroms} 2>{log}"""


rule read_counter:
    input:
        bam="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
        bai="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
    output:
        wig="results/{ID}/icorCNA/readcounts/{SAMPLE}_processed.{GENOME}.wig",
    log:
        "results/logs/{ID}/mapping/{SAMPLE}_all.{GENOME}.log",
    params:
        window=1000000,
        quality=20,
        chroms="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY",
    conda:
        "../envs/icorCNA.yaml"
    shell:
        """
        readCounter --window {params.window} --quality {params.quality} \
        --chromosome {params.chroms} {input.bam} 1>{output.wig} 2>{log}
        """


rule ichorCNA:
    input:
        bam="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
        bai="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
        wig="results/{ID}/icorCNA/readcounts/{SAMPLE}_processed.{GENOME}.wig",
    output:
        outDir=directory("results/{ID}/icorCNA/{SAMPLE}_processed_{GENOME}"),
    params:
        sID="{SAMPLE}_processed",
        ploidy='"c(2,3)"',
        normal='"c(0.5,0.6,0.7,0.8,0.9)"',
        maxCN=5,
        gcWig="../ichorCNA/inst/extdata/gc_hg38_1000kb.wig",
        mapWIG="../ichorCNA/inst/extdata/map_hg38_1000kb.wig",
        centro="../ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt",
        normalPanel="../ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds",
        includeHOMD="False",
        chrs='"c(1:22, "X")"',
        chrTrain='"c(1:22)"',
        estimateNormal="True",
        estimatePloidy="True",
        estimateScPrevalence="True",
        scStates='"c(1,3)"',
        txnE=0.9999,
        txnStrength=10000,
    log:
        "results/logs/{ID}/mapping/{SAMPLE}_all.{GENOME}.log",
    conda:
        "../envs/icorCNA.yaml"
    shell:
        """
        runIchorCNA.R --id {params.sID} --WIG {input.wig} \
        --ploidy {params.ploidy} --normal {params.normal} --maxCN {params.maxCN} \
        --gcWig {params.gcWIG} --mapWig {params.mapWIG} \
        --centromere {params.centro} --normalPanel {params.normalPanel} \
        --includeHOMD {params.includeHOMD} --chrs {params.chrs} \
        --chrTrain {params.chrTrain} --estimateNormal {params.estimateNormal} \
        --estimatePloidy {params.estimatePloidy} \
        --estimateScPrevalence {params.estimateScPrevalence} --scStates {params.scStates} \
        --txnE {params.txnE} --txnStrength {params.txnStrength} --outDir {output.outDir}
        """
