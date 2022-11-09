
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
        window = config["read_counter"]["window"],
        quality = config["read_counter"]["quality"],
        chroms= config["read_counter"]["chroms"],
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
        resources="resources/ichorCNA/"
    output:
        outDir=directory("results/{ID}/icorCNA/{SAMPLE}_processed_{GENOME}/"),
        summary="results/{ID}/icorCNA/{SAMPLE}_processed_{GENOME}/{SAMPLE}_processed.params.txt"
    params:
        sID="{SAMPLE}_processed",
        ploidy = config["ichorCNA"]["ploidy"],
        normal = config["ichorCNA"]["normal"],
        maxCN = config["ichorCNA"]["maxCN"],
        gcWIG = lambda wc: config["ichorCNA"]["gcWIG"][wc.GENOME],
        mapWIG = lambda wc: config["ichorCNA"]["mapWIG"][wc.GENOME],
        centro = lambda wc: config["ichorCNA"]["centro"][wc.GENOME],
        normalPanel = config["ichorCNA"]["normalPanel"],
        includeHOMD = config["ichorCNA"]["includeHOMD"],
        chrs = config["ichorCNA"]["chrs"],
        chrTrain = config["ichorCNA"]["chrTrain"],
        estimateNormal = config["ichorCNA"]["estimateNormal"],
        estimatePloidy = config["ichorCNA"]["estimatePloidy"],
        estimateScPrevalence = config["ichorCNA"]["estimateScPrevalence"],
        scStates = config["ichorCNA"]["scStates"],
        txnE = config["ichorCNA"]["txnE"],
        txnStrength = config["ichorCNA"]["txnStrength"],
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
