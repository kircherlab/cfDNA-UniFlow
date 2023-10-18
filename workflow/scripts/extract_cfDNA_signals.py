#!/usr/bin/env python

import argparse
import gzip
import math
import os
import random
import sys
from collections import defaultdict
import re

import click
import numpy as np
import pandas as pd
import py2bit
import pybedtools as pbt
import pysam
from bx.intervals.intersection import Intersecter, Interval
from cfDNA_GCcorrection import bamHandler
from cfDNA_GCcorrection.utilities import (getGC_content, hash_file, map_chroms,
                                          read_precomputed_table)
from csaps import csaps
from loguru import logger
from matplotlib import pyplot as plt
from mpire import WorkerPool
from scipy import stats
from scipy.signal import savgol_filter
from scipy.ndimage import median_filter

### define functions ###


def isSoftClipped(cigar):
    # Op BAM Description
    # M 0 alignment match (can be a sequence match or mismatch)
    # I 1 insertion to the reference
    # D 2 deletion from the reference
    # N 3 skipped region from the reference
    # S 4 soft clipping (clipped sequences present in SEQ)
    # H 5 hard clipping (clipped sequences NOT present in SEQ)
    # P 6 padding (silent deletion from padded reference)
    # = 7 sequence match
    # X 8 sequence mismatch
    for (op, count) in cigar:
        if op in [4, 5, 6]:
            return True
    return False


def aln_length(cigarlist):
    tlength = 0
    for operation, length in cigarlist:
        if operation == 0 or operation == 2 or operation == 3 or operation >= 6:
            tlength += length
    return tlength


STANDARD_CHROMOSOMES = (
    [str(i) for i in range(1, 23)]
    + ["X", "Y"]
    + ["chr" + str(i) for i in range(1, 23)]
    + ["chrX", "chrY"]
)


def roundGCLenghtBias(gc):
    gc_frac, gc_int = math.modf(round(gc * 100, 2))
    gc_new = gc_int + rng.binomial(1, gc_frac)
    return int(gc_new)


def filter_read(read):
    if read.is_duplicate or read.is_qcfail or read.is_unmapped:
        return None
    if isSoftClipped(read.cigar):
        return None
    if read.is_paired:
        if read.mate_is_unmapped:
            return None
        if read.next_reference_id != read.reference_id:
            return None
    return read


def calculate_signals(
    chrom, start, end, bam, reference, chr_name_bam_to_bit_mapping, R_gc_dict=None, name="unnamed_ROI", strand=".", verbose=False, downsample=None, minInsSize = None, maxInsSize = None, protection=60, gc_correct=False, merged=False, trimmed=False, lengthSR=76 
):
    bam = bamHandler.openBam(bam)
    reference = py2bit.open(reference)
    try:
        tbit_chrom = chr_name_bam_to_bit_mapping[chrom]
    except KeyError as e:
        logger.error(f"Chromosome {chrom} not found in reference mapping. Skipping: {chrom}:{start}-{end}.")
        return
    #logger.debug(f"chrom: {chrom}; mapped tbit_chrom: {tbit_chrom}")
    ID = f"{name}"
    coordinate = f"{chrom}:{start}-{end}"
    regionStart,regionEnd = int(start),int(end)
    tag = 1
    # posRange is a dict of lists of two values. 
    # Each key is a position, the values are coverage, endpoints, average value, number of reads.
    # The dict get updated read per read and position per position.
    posRange = defaultdict(lambda: [0, 0, 0, 0]) #defaultdict(lambda: [0, 0, 0])#
    # filteredReads is an instance of an interval tree (intersecter) class of the bx python package.
    # It is a data structure for performing intersect queries on a set of intervals
    # In this case, intervals represent reads that passed the filters and contain additional information
    filteredReads = Intersecter()
    
    for read in bam.fetch(chrom, max(regionStart - protection - 1, 0), regionEnd + protection + 1,):
        read = filter_read(read)
        if not read:
            continue
        elif read.is_paired:
            if read.is_read1 or (read.is_read2 and read.pnext + read.qlen < regionStart - protection - 1):
                if read.isize == 0:
                    continue
                if downsample != None and random.random() >= downsample:
                    continue
                rstart = min(read.pos, read.pnext) + 1  # 1-based
                lseq = abs(read.isize)
                r_len = abs(read.template_length)
                rend = rstart + lseq - 1  # end included
                if minInsSize != None and ((lseq < minInsSize) or (lseq > maxInsSize)):
                    continue
            else:
                continue
        else:
            if downsample != None and random.random() >= downsample:
                continue
            rstart = read.pos + 1  # 1-based
            lseq = aln_length(read.cigar)
            r_len = read.query_length
            rend = rstart + lseq - 1  # end included
            #logger.debug(rstart, rend)
            if minInsSize != None and ((lseq < minInsSize) or (lseq > maxInsSize)):
                continue
        # if GC should be used, change the value of the tag to the corresponding bias value
        
        if gc_correct:
            try:
                gc = getGC_content(
                    reference,
                    tbit_chrom,
                    read.reference_start,
                    read.reference_end,
                    fraction=True,
                )
                gc = roundGCLenghtBias(gc)
            except Exception as detail:
                if verbose:
                    logger.exception(detail)
                continue
            if R_gc_dict:
                try:
                    tag = round(float(1) / R_gc_dict[r_len][str(gc)],5)#2) ; should be set by round_digit
                except KeyError as e:
                    tag = 1
            else:
                try:
                    tag=round(dict(read.tags)["YC"],5)
                except KeyError as e:
                    tag=1
        # add tag to interval. tag defaults to 1 if there is no YC weight attached to the read or gc_correct is falsy
        filteredReads.add_interval(Interval(rstart, rend, value={"YC": tag}))
        if read.is_paired:
            for i in range(rstart, rend + 1):
                if i >= regionStart and i <= regionEnd:
                    posRange[i][0] += tag
                    posRange[i][2] += 1
                    posRange[i][3] += (tag-posRange[i][3])/posRange[i][2]
                if rstart >= regionStart and rstart <= regionEnd:
                    posRange[rstart][1] += tag  
                if rend >= regionStart and rend <= regionEnd:
                    posRange[rend][1] += tag
        else:
            for i in range(rstart, rend + 1):
                if i >= regionStart and i <= regionEnd:
                    posRange[i][0] += tag
                    posRange[i][2] += 1
                    posRange[i][3] += (tag-posRange[i][3])/posRange[i][2]
            if (merged or read.qname.startswith("M_")) or (
                (trimmed or read.qname.startswith("T_"))
                and read.qlen <= lengthSR - 10
            ):
                if rstart >= regionStart and rstart <= regionEnd:
                    posRange[rstart][1] += tag  
                if rend >= regionStart and rend <= regionEnd:
                    posRange[rend][1] += tag  
            elif read.is_reverse:
                if rend >= regionStart and rend <= regionEnd:
                    posRange[rend][1] += tag  
            else:
                if rstart >= regionStart and rstart <= regionEnd:
                    posRange[rstart][1] += tag  
        
    logger.debug("Evaluating posRange vector...\n")
    # filename = outfile%cid
    # outfile = gzip.open(filename,'w')
    cov_sites = 0
    outLines = []
    wps_list = []
    cov_list = []
    starts_list = []
    mean_weight_list = []
    
    for pos in range(regionStart, regionEnd + 1):
        rstart, rend = pos - protection, pos + protection
        gcount, bcount, ecount = 0.0, 0.0, 0.0
        for read in filteredReads.find(rstart, rend):
                if (read.start > rstart) or (read.end < rend):
                    bcount += read.value["YC"]
                else:
                    gcount += read.value["YC"]

        covCount, startCount, n_reads, mean_value = posRange[pos]
        #logger.debug(posRange[pos])
        #covCount, startCount, n_reads = posRange[pos]
        #mean_value = covCount/n_reads
        cov_sites += covCount
        wpsValue = gcount - (bcount + ecount)
        
        
        wps_list.append(wpsValue)
        cov_list.append(covCount)
        starts_list.append(startCount)
        mean_weight_list.append(mean_value)
    
    if strand == "-":
        wps_list = wps_list[::-1]
        cov_list = cov_list[::-1]
        starts_list = starts_list[::-1]
        mean_weight_list = mean_weight_list[::-1]

    
    return ID,coordinate,wps_list,cov_list,starts_list,mean_weight_list


@click.command()
@click.option(
    "--bamfile",
    "-b",
    "bam_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Sorted BAM file.",
)
@click.option(
    "--genome",
    "-g",
    "reference_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="""Genome reference in two bit format. Most genomes can be 
            found here: http://hgdownload.cse.ucsc.edu/gbdb/ 
            Search for the .2bit ending. Otherwise, fasta 
            files can be converted to 2bit using the UCSC 
            programm called faToTwoBit available for different 
            plattforms at '
            http://hgdownload.cse.ucsc.edu/admin/exe/""",
)
@click.option(
    "--regions",
    "-r",
    "regions",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Regions of interest in .bed fromat.",
)
@click.option(
    "--outfile",
    "-o",
    "outfile",
    show_default=True,
    type=click.Path(writable=True),
    help="""Outfile prefix to save the extracted signals.
            The filename must containt %s or {}, which will be substituted with "WPS" and "COV".
            Output file will be gzipped tab separated files.""",
    default="signal_{}.tsv.gz",
)
@click.option(
    "-gc",
    "--gc_correct",
    "GC_correct",
    is_flag=True,
    help="""Flag: indacte if gc correction should be used. 
            Expects bam file with GC weights attached to reads as YC tags. 
            Is automatically set to true if a GCbiasFrequencyfile is supplied.""",
)
@click.option(
    "--GCbiasFrequenciesFile",
    "-freq",
    "GC_file",
    type=click.Path(exists=True, readable=True),
    help="""Indicate the output file from computeGCBias containing 
            the coverage bias values per GC-content and fragment length.
            This file should be a (gzipped) tab separated file.""",
)
@click.option(
    "--minInsSize",
    "-min",
    "minInsSize",
    default=None,
    type=click.INT,
    help="Minimum read length threshold to consider for signal extraction.",
)
@click.option(
    "--maxInsSize",
    "-max",
    "maxInsSize",
    default=None,
    type=click.INT,
    help="Maximum read length threshold to consider for signal extraction.",
)
@click.option(
    "--max_length",
    "max_length",
    default=1000,
    type=click.INT,
    show_default=True,
    help="Assumed maximum insert size.",
)
@click.option(
    "--protection",
    "-w",
    "protection",
    default=120,
    type=click.INT,
    help="Base pair protection window size.",
)
@click.option(
    "--lengthSR",
    "-l",
    "lengthSR",
    default=76,
    type=click.INT,
    help="Length of full reads.",
)
@click.option(
    "--merged",
    "-m",
    "merged",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Assume reads are merged.""",
)
@click.option(
    "--trimmed",
    "-t",
    "trimmed",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Assume reads are trimmed.""",
)
@click.option(
    "--num_cpus",
    "-p",
    "num_cpus",
    type=click.INT,
    default=1,
    show_default=True,
    help="Number of processors to use.",
)
@click.option(
    "--seed",
    "seed",
    default=None,
    type=click.INT,
    help="""Set seed for reproducibility.""",
)
@click.option("-v", "--verbose", "verbose", is_flag=True, help="Enables verbose mode")
@click.option("--debug", "debug", is_flag=True, help="Enables debug mode")
@click.option(
    "--filter_size",
    "filter_size",
    default=None,
    type=click.INT,
    help="Filter size.",
)
@click.option(
    "--round",
    "round_digits",
    default=2,
    type=click.INT,
    help="Digits to round the correction values to.",
)
@click.option(
    "-sc",
    "--standard_chroms",
    "standard_chroms",
    is_flag=True,
    help="Flag: filter chromosomes to human standard chromosomes.",
)
def main(
    bam_file,
    reference_file,
    regions,
    outfile,
    GC_correct,
    GC_file,
    minInsSize,
    maxInsSize,
    max_length,
    protection,
    lengthSR,
    merged,
    trimmed,
    num_cpus,
    seed,
    verbose,
    debug,
    filter_size,
    round_digits,
    standard_chroms
):

    ### initial setup
    progress_bar = False

    if debug:
        debug_format = "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> | <level>process: {process}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"
        logger.remove()
        logger.add(
            sys.stderr, level="DEBUG", format=debug_format, colorize=True, enqueue=True
        )
        logger.debug("Debug mode active.")
    elif verbose:
        progress_bar = True
        info_format = (
            "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{message}</level>"
        )
        logger.remove()
        logger.add(
            sys.stderr, level="INFO", format=info_format, colorize=False, enqueue=True
        )
    else:
        logger.remove()
        logger.add(
            sys.stderr,
            level="WARNING",
            colorize=True,
            enqueue=True,
        )

    global rng
    rng = np.random.default_rng(seed=seed)
    random.seed(seed)

    ### execute Utility functions

    outfile = outfile.replace("%s", "{}")

    tbit = py2bit.open(reference_file)

    bam, mapped, unmapped, stats = bamHandler.openBam(
        bam_file, returnStats=True, nThreads=num_cpus
    )

    chr_name_bam_to_bit_mapping = map_chroms(
        bam.references,
        list(tbit.chroms().keys()),
        ref_name="bam file",
        target_name="2bit reference file",
    )

    protection = protection // 2

    ### parse regions in .bed format and create an iterator
    regions = pbt.BedTool(regions).to_dataframe(
            dtype={"chrom": "category", "start": "uint32", "end": "uint32"}
        ).drop(columns=["score"])
    if standard_chroms:
        regions = regions.loc[regions["chrom"].isin(STANDARD_CHROMOSOMES)]
    regions_to_ref_mapping = map_chroms(set(regions["chrom"].unique()), list(tbit.chroms().keys()))
    regions = regions.replace({"chrom":regions_to_ref_mapping})
    region_iterator = regions.to_dict(orient="records")


    if GC_file:
        data = pd.read_csv(GC_file, sep="\t", index_col=[0, 1])
        R_gc = data.loc["R_gc"]
        if filter_size:
            logger.info(f"Using filtersize {filter_size} for GC bias smoothing.")
            tmp_index = R_gc.index
            tmp_cols = R_gc.columns
            R_gc = pd.DataFrame(median_filter(R_gc,size=filter_size))
            R_gc.index = tmp_index
            R_gc.columns = tmp_cols

        logger.info(R_gc.T.describe().T.describe())

        R_gc_dict = R_gc.round(5).to_dict(orient="index")

    ### setup parameter dictionary ###

    param_dict = {
        "bam": bam_file,
        "reference": reference_file,
        "chr_name_bam_to_bit_mapping": chr_name_bam_to_bit_mapping,
        "minInsSize": minInsSize,
        "maxInsSize": maxInsSize,
    }

    ### add values to dict, if they deviate from defaults
    if GC_file:
        param_dict["gc_correct"] = True
        param_dict["R_gc_dict"] = R_gc_dict
    if GC_correct:
        param_dict["gc_correct"] = True
    if verbose or debug:
        param_dict["verbose"] = True
    if protection != 60:
        param_dict["protection"] = protection
    if lengthSR != 76:
        param_dict["lengthSR"] = lengthSR
    if merged:
        param_dict["merged"] = True
    if trimmed:
        param_dict["trimmed"] = True
    if max_length != 1000:
        param_dict["max_length"] = max_length

    tasks = ({**region, **param_dict} for region in region_iterator)

    ### setup worker pool and distribute tasks

    with WorkerPool(n_jobs=num_cpus) as pool:
        imap_res = pool.imap(calculate_signals, tasks, iterable_len=len(region_iterator), progress_bar=progress_bar)

    ### stream data to output files ###

    # It would be nice to have an open function that automatically determines if file should be gzipped
    counter = 0
    with gzip.open(outfile.format("WPS"), "wt") as wps, gzip.open(
        outfile.format("COV"), "wt"
    ) as cov, gzip.open(
        outfile.format("MEAN_WEIGHT"), "wt"
    ) as mean_weight:
        for res in imap_res:
            if not res:
                continue
            else:
                ID,coordinate, wps_list, cov_list, starts_list, mean_weight_list = res
            if counter % 1000 == 0:
                logger.debug(f"writing data for: {ID} {coordinate}")
            wps.write(
                ID
                + ","
                + coordinate
                + ","
                + ",".join(map(lambda x: re.sub('.0$','',str(round(x, 5))), wps_list))
                + "\n"
            )
            cov.write(
                ID
                + ","
                + coordinate
                + ","
                + ",".join(map(lambda x: re.sub('.0$','',str(round(x, 5))), cov_list))
                + "\n"
            )
            mean_weight.write(
                ID
                + ","
                + coordinate
                + ","
                + ",".join(map(lambda x: re.sub('.0$','',str(round(x, 5))), mean_weight_list))
                + "\n"
            )
            counter += 1


if __name__ == "__main__":
    main()
