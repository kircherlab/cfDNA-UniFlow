def get_chroms_readcounter(chromfile):
    if config["read_counter"]["chroms"]["autoselect_chroms"]:
        with open(chromfile, "r") as f:
            chrom = f.read()
        return chrom
    else:
        return config["read_counter"]["chroms"]["chrom_list"]


rule get_chroms:
    input:
        bam="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
        bai="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
    output:
        chroms="results/{ID}/icorCNA/chroms/{SAMPLE}_processed.{GENOME}.chromosomes.txt",
    params:
        strict_filter="awk '{if ($3 > 0) {printf $0 \"\\n\"}}' | "
        if config["read_counter"]["chroms"]["autoselect_strict"]
        else "",
    log:
        "logs/{ID}/get_chroms/get_chroms_{SAMPLE}_processed.{GENOME}.chromosomes.log",
    conda:
        "../envs/cfDNA_prep.yaml"
    shell:
        r"""(samtools idxstats {input.bam} | \
        {params.strict_filter} \
        cut -f 1 | grep -w 'chr[1-9]\|chr[1-2][0-9]\|chr[X,Y]\|^[1-2][0-9]\|^[0-9]\|^[X,Y]' \
        | sort -k1,1 -V -s | tr '\n' ','| sed 's/,*\r*$//' \
        1> {output.chroms}) 2>{log}"""


rule read_counter:
    input:
        bam="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam",
        bai="results/{ID}/mapped_reads/{SAMPLE}_processed.{GENOME}.bam.bai",
        chroms="results/{ID}/icorCNA/chroms/{SAMPLE}_processed.{GENOME}.chromosomes.txt",
    output:
        wig="results/{ID}/icorCNA/readcounts/{SAMPLE}_processed.{GENOME}.wig",
    log:
        "logs/{ID}/read_counter/read_counter_{SAMPLE}_processed.{GENOME}.log",
    params:
        window=config["read_counter"]["window"],
        quality=config["read_counter"]["quality"],
        chroms=lambda wc: get_chroms_readcounter(
            f"results/{wc.ID}/icorCNA/chroms/{wc.SAMPLE}_processed.{wc.GENOME}.chromosomes.txt"
        ),
    conda:
        "../envs/ichorCNA.yaml"
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
        resources="resources/ichorCNA/",
        chroms="results/{ID}/icorCNA/chroms/{SAMPLE}_processed.{GENOME}.chromosomes.txt",
    output:
        outDir=directory("results/{ID}/icorCNA/{SAMPLE}_processed_{GENOME}/"),
        plotDir=report(
            directory(
                "results/{ID}/icorCNA/{SAMPLE}_processed_{GENOME}/{SAMPLE}_processed/"
            ),
            patterns=[
                "{SAMPLE}_processed_genomeWide.png",
            ],
            caption="../report/ichorCNA.rst",
            category="ichorCNA",
            labels={"Sample": "{SAMPLE}", "Type": "CNA-TumorFraction plot"},
        ),
        summary="results/{ID}/icorCNA/{SAMPLE}_processed_{GENOME}/{SAMPLE}_processed.params.txt",
    log:
        "logs/{ID}/ichorCNA/ichorCNA_{SAMPLE}_processed.{GENOME}.log",
    params:
        genome="{GENOME}",
        plot_format="png",
        sID="{SAMPLE}_processed",
        ploidy=config["ichorCNA"]["ploidy"],
        normal=config["ichorCNA"]["normal"],
        maxCN=config["ichorCNA"]["maxCN"],
        gcWIG=lambda wc: config["ichorCNA"]["gcWIG"][wc.GENOME],
        mapWIG=lambda wc: config["ichorCNA"]["mapWIG"][wc.GENOME],
        centro=lambda wc: config["ichorCNA"]["centro"][wc.GENOME],
        normalPanel=lambda wc: config["ichorCNA"]["normalPanel"][wc.GENOME],
        includeHOMD=config["ichorCNA"]["includeHOMD"],
        chrs=config["ichorCNA"]["chrs"],
        chrTrain=config["ichorCNA"]["chrTrain"],
        estimateNormal=config["ichorCNA"]["estimateNormal"],
        estimatePloidy=config["ichorCNA"]["estimatePloidy"],
        estimateScPrevalence=config["ichorCNA"]["estimateScPrevalence"],
        scStates=config["ichorCNA"]["scStates"],
        txnE=config["ichorCNA"]["txnE"],
        txnStrength=config["ichorCNA"]["txnStrength"],
    conda:
        "../envs/ichorCNA.yaml"
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
        --txnE {params.txnE} --txnStrength {params.txnStrength} --plotFileType {params.plot_format} \
        --outDir {output.outDir} 2>{log}
        """
