# cfDNA UniFlow - Unified Preprocessing Pipeline for cell-free DNA from liquid biopsies

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.4.1-brightgreen.svg)](https://snakemake.bitbucket.io)
<div align="justify">

cfDNA UniFlow is a unified, standardized, and ready-to-use workflow for processing WGS cfDNA samples from liquid biopsies. It includes essential steps for pre-processing raw cfDNA samples, quality control and reporting. Additionally, several optional utility functions like GC bias correction and estimation of copy number state are included. Figure S1 gives a detailed overview of the workflow.
</div>

<figure>
 <img loading="lazy" src="supplement/cfDNA_unifyed_preprocessing.drawio.png">
 <figcaption>
 <div align="justify">
  <strong>Figure S1: Overview of cfDNA Uniflow.</strong>  Functionalities are color coded by task. Red boxes represents rules for the automatic download of public resources. Grey boxes are optional steps. Blue boxes containt the core functionailty of cfDNA Uniflow. Green boxes are optional, but highly recommended steps and yellow boxes summarize the Quality Control and reporting steps.
 </div>
 </figcaption>
</figure>

## Authors <!-- omit from toc -->

* Sebastian Röner (@sroener)

## Table of Contents <!-- omit from toc -->

- [cfDNA UniFlow - Unified Preprocessing Pipeline for cell-free DNA from liquid biopsies](#cfdna-uniflow---unified-preprocessing-pipeline-for-cell-free-dna-from-liquid-biopsies)
    - [1 Dependencies](#1-dependencies)
    - [2 Setup](#2-setup)
        - [Step 1: Obtain a copy of this workflow](#step-1-obtain-a-copy-of-this-workflow)
        - [Step 2: Install Snakemake](#step-2-install-snakemake)
        - [Step 3: Configure Workflow](#step-3-configure-workflow)
        - [Step 4: Execute workflow](#step-4-execute-workflow)
        - [Step 5: Investigate results](#step-5-investigate-results)
    - [3 Functional summary](#3-functional-summary)
        - [3.1 Raw data processing](#31-raw-data-processing)
        - [3.2 Quality control](#32-quality-control)
        - [3.3 Utility functionality](#33-utility-functionality)
        - [3.4 Report](#34-report)
        - [3.5 Notes on resource requirements](#35-notes-on-resource-requirements)
    - [4 Quickstart guide](#4-quickstart-guide)
        - [4.1 Download test files](#41-download-test-files)
        - [4.2 Check config files](#42-check-config-files)
        - [4.3 Executing the workflow](#43-executing-the-workflow)

<div align="justify">

## 1 Dependencies

To minimize conflicts between packages, all dependencies for the workflow rules are managed via separate conda environments located in `./workflow/envs`. They get automatically installed and used when Snakemake is executed with the `--use-conda` flag. This is the recommended way of executing the workflow.

The only exception is NGmerge, a read merging and adapter removal program, which is included in the GitHub repository due to an outdated bioconda recipe. The NGmerge executable was downloaded and compiled as described in the official documentation of NGmerge v0.3 and is located in the scripts directory (./workflow/scripts/NGmerge). Additionally, we provide an adjusted quality profile for Phred+33 scores of Illumina 1.8+, which ranges from 0 to 41 instead of 0 to 40 in earlier versions. The quality profile file was modified by duplicating the last column and appending it as a new column.

## 2 Setup

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Install Snakemake

For best compatibility, it is recommended to execute cfDNA UniFlow via conda and Snakemake. For this, it is required to first install conda. If conda is not installed yet, follow the [official documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html).

After successful installation, set up an environment for Snakemake. This can be done by executing the following command:

```bash
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

The environment can be activated via the `conda activate snakemake` command.

For installation details, see the instructions in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 3: Configure Workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` according to the included comments to configure the workflow execution, and `samples.tsv` to specify your sample setup. Both files will be automatically validated before workflow execution.

### Step 4: Execute workflow

Activate the conda environment:

```bash
conda activate snakemake
```

Test your configuration by performing a dry-run with `<CONFIGFILE>` representing the path to your modified config file via:

```bash
snakemake --use-conda --configfile <CONFIGFILE> -n
```

Execute the workflow locally using `$N` cores:

```bash
snakemake --use-conda --configfile <CONFIGFILE> --cores $N
```

See the Snakemake documentation on [workflow execution](https://snakemake.readthedocs.io/en/stable/executing/cli.html) and execution in [cluster environments](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) for further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

```bash
snakemake --configfile <CONFIGFILE> --report report.html
```

This report can, e.g., be forwarded to your collaborators.

A functional description of reporting can be found in [section 3.4](#34-report).

An example report in zip format can be found in the supplement directory. For viewing the zip file needs to be extracted.

## 3 Functional summary

### 3.1 Raw data processing

The core functionality of cfDNA UniFlow is the processing of Whole Genome Sequencing (WGS) data from liquid biopsies. Input data is expected in either FASTQ or BAM format. If a BAM file was provided, it gets converted to FASTQ files using SAMTools. Afterwards, several steps for improving read quality and preparation for mapping are executed, for which two options are provided. Either the recommended merging of reads with NGmerge, which removes adapters and corrects sequencing errors, or trimming of adapters with Trimmomatic. In both cases, results are filtered for a specified minimum read length. Remaining reads are mapped to a reference genome via bwa-mem2. If NGmerge was used for adapter removal, it is possible to include reads in mapping for which only adapters were removed when merging was not possible. The core processing is finalized by marking duplicates and creating a bam index with SAMTools markdup and index. Processed reads are then submitted for Quality Control and optional steps.

### 3.2 Quality control

In the quality control step, general post-alignment statistics and graphs are calculated for each sample via SAMTools stats and FastQC. Additionally, sample-wide coverage statistics and coverage at different genomic regions are calculated with Mosdepth, a fast BAM/CRAM depth calculation tool for WGS, exome, or targeted sequencing. The QC results are aggregated in HTML report via MultiQC.

### 3.3 Utility functionality

In addition to the preprocessing and quality control functionality, cfDNA UniFlow contains some utility functions. The first is the widely used tool IchorCNA, which can be used for predicting copy number alteration (CNA) states across the genome. Further, it uses this information for estimating tumor fractions in cfDNA samples.

The second utility function is our [inhouse GC bias estimation method](https://github.com/kircherlab/cfDNA_GCcorrection). It can not only be used for estimating fragment length and GC-content dependent technical biases, but also includes the option of attaching correction values to the reads. These can be used downstream for a wide variety of signal extraction methods, while preserving the original read coverage patterns.

### 3.4 Report

Finally, all results and summary statistics for the specified samples are aggregated in one report, making a wide variety of information easily accessible. The report file is generated using Snakemake’s report feature. After the workflow finished, the report can be generated by executing `Snakemake –configfile <CONFIGFILE> --snakefile <SNAKEFILE> --report <REPORTNAME>.html`. For better report performance, it is recommended to use `--report <REPORTNAME>.zip`, which creates a zipped directory structure with the needed information instead of saving it in the HTML file itself. More information on reporting can be found in the official [Snakemake  documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html).

### 3.5 Notes on resource requirements

* The index creation of bwa-mem2 is resource intensive:

```bash
# Indexing the reference sequence (Requires 28N GB memory where N is the size of the reference sequence).
./bwa-mem2 index [-p prefix] <in.fasta>
Where 
<in.fasta> is the path to reference sequence fasta file and 
<prefix> is the prefix of the names of the files that store the resultant index. Default is in.fasta.
```

* bwa-mem2 mem uses around 4GB memory per thread

## 4 Quickstart guide

Our goal in developing cfDNA UniFlow is to provide a scalable, configurable and easy-to-use workflow specifically tailored towards the processing of cfDNA samples. Users only need to provide sequencing information in FASTQ or BAM format and optionally modify the configuration file to their needs. Here we provide an example with a small input file for testing the workflows functionality.

**Note:** For simplicity, we expect a Unix system for this guide.

### 4.1 Download test files

First, you should create the expected directory for the test files:

```bash
mkdir -p resources/testsample/
```

Afterwards, download the files located on this [webserver](https://kircherlab.bihealth.org/download/cfDNA-testSam):

```bash
curl -L -o resources/testsample/testsample_hg19_1x_chr20-22.bam https://kircherlab.bihealth.org/download/cfDNA-testSample/testsample_hg19_1x_chr20-22.bam
```

### 4.2 Check config files

There are two files that are used in this example:

- config/test-config.yaml
- config/test-samples.tsv

Both files don't need to be edited and are configured for a quick functionality test.

### 4.3 Executing the workflow

The last step is executing the Workflow. For this you need to be in the root directory of the cloned repository.

Test the configuration with a dry-run:

```bash
snakemake --use-conda --configfile config/test-config.yaml  -n
```

The workflow is executed locally with $N cores via:

```bash
snakemake --use-conda --configfile config/test-config.yaml --cores $N
```

For cluster execution, read the guidelines in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).


</div>
