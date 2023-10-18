# cfDNA UniFlow: A unified preprocessing pipeline for cell-free DNA from liquid biopsies
<div align="justify">

cfDNA UniFlow is a unified, standardized, and ready-to-use workflow for processing whole genome sequencing (WGS) cfDNA samples from liquid biopsies. It includes essential steps for pre-processing raw cfDNA samples, quality control and reporting. Additionally, several optional utility functions like GC bias correction and estimation of copy number state are included. Finally, we provide specialized methods for extracting coverage derived signals and visualizations comparing cases and controls. Figure S1 gives a detailed overview of the workflow.
</div>

<figure>
 <img loading="lazy" src="supplement/cfDNA_unifyed_preprocessing.drawio.png">
 <figcaption>
 <div align="justify">
  <strong>Figure S1: Overview of cfDNA Uniflow.</strong> Functionalities are color coded by task. Red boxes represent rules for the automatic download of public resources. Grey boxes are optional steps. Blue boxes contain the core functionalty of cfDNA Uniflow. Green boxes are optional, but highly recommended steps and yellow boxes summarize the Quality Control and reporting steps.
 </div>
 </figcaption>
</figure>

## Authors <!-- omit from toc -->

* Sebastian Röner (@sroener)

<div style="page-break-after: always;"></div>

## Table of Contents <!-- omit from toc -->

- [cfDNA UniFlow: A unified preprocessing pipeline for cell-free DNA from liquid biopsies](#cfdna-uniflow-a-unified-preprocessing-pipeline-for-cell-free-dna-from-liquid-biopsies)
    - [1 Dependencies](#1-dependencies)
    - [2 Functional summary](#2-functional-summary)
        - [2.1 Raw data processing](#21-raw-data-processing)
        - [2.2 Quality control](#22-quality-control)
        - [2.3 Signal extraction and utility functionality](#23-signal-extraction-and-utility-functionality)
        - [2.4 Report](#24-report)
        - [2.5 Notes on resource requirements](#25-notes-on-resource-requirements)
    - [3 Quick-start Guide](#3-quick-start-guide)
        - [3.1 Setup](#31-setup)
            - [Step 1: Obtain a copy of this workflow](#step-1-obtain-a-copy-of-this-workflow)
            - [Step 2: Install Snakemake](#step-2-install-snakemake)
            - [Step 3: Activate conda environment](#step-3-activate-conda-environment)
        - [3.2 Integration test](#32-integration-test)
            - [Step 1: Download test files](#step-1-download-test-files)
            - [Step 2: Check config files](#step-2-check-config-files)
            - [Step 3: Executing the workflow](#step-3-executing-the-workflow)
        - [3.3 Running the workflow with your own data](#33-running-the-workflow-with-your-own-data)
            - [Step 1: Configure Workflow](#step-1-configure-workflow)
            - [Step 2: Execute workflow](#step-2-execute-workflow)
            - [Step 3: Investigate results](#step-3-investigate-results)

<div style="page-break-after: always;"></div>

<div align="justify">

## 1 Dependencies

To minimize conflicts between packages, all dependencies for the workflow rules are managed via separate conda environments located in `./workflow/envs`. They get automatically installed and used when Snakemake is executed with the `--use-conda` flag. This is the recommended way of executing the workflow.

The only exception is NGmerge, a read merging and adapter removal program, which is included in the GitHub repository due to an outdated bioconda recipe. The NGmerge executable was downloaded and compiled as described in the official documentation of NGmerge v0.3 and is located in the scripts directory (./workflow/scripts/NGmerge). Additionally, we provide an adjusted quality profile for Phred+33 scores of Illumina 1.8+, which ranges from 0 to 41 instead of 0 to 40 in earlier versions. The quality profile file was modified by duplicating the last column and appending it as a new column.

## 2 Functional summary

### 2.1 Raw data processing

The core functionality of cfDNA UniFlow is the processing of Whole Genome Sequencing (WGS) data from liquid biopsies. Input data is expected in either FASTQ or BAM format. If a BAM file was provided, it gets converted to FASTQ files using SAMtools. Afterwards, several steps for improving read quality and preparation for mapping are executed, for which two options are provided. Either the recommended merging of reads with NGmerge, which removes adapters and corrects sequencing errors, or trimming of adapters with Trimmomatic. In both cases, results are filtered for a specified minimum read length. Remaining reads are mapped to a reference genome via bwa-mem2. If NGmerge was used for adapter removal, it is possible to include reads in mapping for which only adapters were removed when merging was not possible due to no sufficient overhang between read pairs. The core processing is finalized by marking duplicates and creating a bam index with SAMtools markdup and index. Processed reads are then submitted for Quality Control and optional steps.

### 2.2 Quality control

In the quality control step, general post-alignment statistics and graphs are calculated for each sample via SAMtools stats and FastQC. Additionally, sample-wide coverage statistics and coverage at different genomic regions are calculated with Mosdepth, a fast BAM/CRAM depth calculation tool for WGS, exome, or targeted sequencing. The QC results are aggregated in HTML report via MultiQC.

<figure>
 <img loading="lazy" src="supplement/example_QC-report.png">
 <figcaption>
 <div align="justify">
  <strong>Figure S2: Example section of a QC report.</strong> The QC report contains post-alignment statistics from SAMtools stats, FastQC and Mosdepth.
 </div>
 </figcaption>
</figure>

### 2.3 Signal extraction and utility functionality

In addition to the preprocessing and quality control functionality, cfDNA UniFlow contains some utility functions. The first is the widely used tool ichorCNA, which can be used for predicting copy number alteration (CNA) states across the genome. Further, it uses this information for estimating tumor fractions in cfDNA samples. By default, we use the recommended settings for ichorCNA, including profiles provided in the ichorCNA repository. However, it is possible to specify custom profiles and parameters in the configuration file.

<figure>
 <img loading="lazy" src="supplement/example_ichorCNA.png">
 <figcaption>
 <div align="justify">
  <strong>Figure S3: ichorCNA plot of copy number alterations.</strong> The plot shows genome wide annotation of estimated copy number alterations, tumor fraction and other parameters used by ichorCNA.This example was generated using a breast cancer sample with an average GC content of 38%.
 </div>
 </figcaption>
</figure>

<div style="page-break-after: always;"></div>

The second utility function is our [inhouse GC bias estimation method](https://github.com/kircherlab/cfDNA_GCcorrection). It can not only be used for estimating fragment length and GC-content dependent technical biases, but also includes the option of attaching correction values to the reads. These can be used downstream for a wide variety of signal extraction methods, while preserving the original read coverage patterns. Additionally, it is possible to include a visualization of the estimated biases in the report.

<figure>
 <img loading="lazy" src="supplement/example_GC-bias.png">
 <figcaption>
 <div align="justify">
  <strong>Figure S4: Global GC bias estimate.</strong> The global GC bias was estimated using ratio between observed and expected fragments stratified by GC content [%] and fragment length [bp].This example was generated using a breast cancer sample with an average GC content of 38%.
 </div>
 </figcaption>
</figure>

Building on the GC bias estimation, we provide a method for extracting coverage derived signals around user defined regions. The resulting signals can be visualized for comparing biased vs corrected states and for comparing cases against controls.

<figure>
 <img loading="lazy" src="supplement/example_GC-overlay.png">
 <figcaption>
 <div align="justify">
  <strong>Figure S5: Regional effects of GC bias correction.</strong> Composite coverage signal of 10000 LYL1 binding sites before and after GC bias correction. Accessibility of the binding sites, indicated by lower coveragem would have been overestimated before correction.  This example was generated using a breast cancer sample with an average GC content of 38%.
 </div>
 </figcaption>
</figure>

<figure>
 <img loading="lazy" src="supplement/example_case-control-GRHL2.png">
 <figcaption>
 <div align="justify">
  <strong>Figure S6: Case-control plot.</strong> GC corrected composite coverage signal of 10000 GRHL2 binding sites in a healthy control and a breast cancer sample. Lower central coverage in the cancer sample is consistent with increased GRHL2 activity in many cancer types.
  </div>
 </figcaption>
</figure>

### 2.4 Report

Finally, all results and summary statistics for the specified samples are aggregated in one report, making a wide variety of information easily accessible. The report file is generated using Snakemake’s report feature. After the workflow finished, the report can be generated by executing `Snakemake –configfile <CONFIGFILE> --snakefile <SNAKEFILE> --report <REPORTNAME>.html`. For better report performance, it is recommended to use `--report <REPORTNAME>.zip`, which creates a zipped directory structure with the needed information instead of saving it in the HTML file itself. More information on reporting can be found in the official [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html).

**Note:** An example report can be found [here](https://github.com/kircherlab/cfDNA-UniFlow/tree/main/supplement)

### 2.5 Notes on resource requirements

* The index creation of bwa-mem2 is resource intensive:

```bash
# Indexing the reference sequence (Requires 28N GB memory where N is the size of the reference sequence).
./bwa-mem2 index [-p prefix] <in.fasta>
Where 
<in.fasta> is the path to reference sequence FASTA file and 
<prefix> is the prefix of the names of the files that store the resultant index. Default is in.FASTA.
```

More information can be found in the [documentation](https://github.com/bwa-mem2/bwa-mem2#usage) or [this issue](https://github.com/bwa-mem2/bwa-mem2/issues/111).

* bwa-mem2 mem uses around 4GB memory per thread.

## 3 Quick-start Guide

### 3.1 Setup

#### Step 1: Obtain a copy of this workflow

1. Create a new GitHub repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 2: Install Snakemake

For best compatibility, it is recommended to execute cfDNA UniFlow via conda and Snakemake. For this, it is required to first install conda. If conda is not installed yet, follow the [official documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html).

After successful installation, set up an environment for Snakemake. This can be done by executing the following command:

```bash
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

The environment can be activated via the `conda activate snakemake` command.

For installation details, see the instructions in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

#### Step 3: Activate conda environment

Activate the conda environment:

```bash
conda activate snakemake
```

### 3.2 Integration test

Our goal in developing cfDNA UniFlow is to provide a scalable, configurable, and easy-to-use workflow specifically tailored towards the processing of cfDNA samples. Users only need to provide sequencing information in FASTQ or BAM format and optionally modify the configuration file to their needs. Here we provide an example with a small input file for testing the workflows functionality.

**Note:** As this is an integration test, the generated data is not intended to generate meaningful results. This test is only meant to check if the workflow is running without errors.

#### Step 1: Download test files

First, you should create the expected directory for the test files:

```bash
mkdir -p resources/testsample/
```

Afterwards, download the files located on this [webserver](https://kircherlab.bihealth.org/download/cfDNA-testSam):

```bash
curl -L -o resources/testsample/testsample_hg19_1x_chr20-22.bam https://kircherlab.bihealth.org/download/cfDNA-testSample/testsample_hg19_1x_chr20-22.bam
```

#### Step 2: Check config files

There are three files that are used in this example:

- config/test-config.yaml
- config/test-samples.tsv
- config/test-regions.tsv

**Note**: These files are configured for a quick integration test and shouldn't be edited.

#### Step 3: Executing the workflow

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

### 3.3 Running the workflow with your own data

After running the integration test, you can start using the workflow with your own data. For this, you need to modify the configuration files and provide your own data.

#### Step 1: Configure Workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` according to the included comments to configure the workflow execution. The most important options are the path to the respective `samples.tsv` and `regions.tsv` files. All other default files are either included, or will be downloaded and generated automatically. 

The `samples.tsv` file should contain the following columns:

```yaml
# Identifier of the cohort/experiment
ID: string
# Sample name/identifier
sample: string
# Path to the BAM file (for reprocessing); insert "-" if FASTQ files are used
bam: string 
# Path to the R1 FASTQ file; insert "-" if bam is used
fq1: string
# Path to the R2 FASTQ file; insert "-" if bam is used
fq2: string
# Name of the target genome build (hg19 or hg38)
genome_build: string
# Name of the library prep (e.g. ThruPLEX DNA-seq)
library_name: string
# Platform used for sequencing (e.g. Illumina NextSeq 500)
platform: string
# Status of the sample (e.g. Case/Control or Healthy/Cancer)
status: string
# Additional information that should be included in read group ID (e.g. condition)
info: string
```

If FASTQ files are used, either single-end or paired-end sequencing can be specified, depending on how many FASTQ files are provided. If BAM files are used, the FASTQ columns should be filled with "-" and vice versa. 

The `regions.tsv` file should contain the following columns:

```yaml
# region name/identifier
target: string
# path to the bed files for regions of interest
path: string
```

**Note:** All three files will be automatically validated before workflow execution.

#### Step 2: Execute workflow

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

#### Step 3: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

```bash
snakemake --configfile <CONFIGFILE> --report report.html
```

This report can, e.g., be forwarded to your collaborators.

A functional description of reporting can be found in [section 2.4](#24-report).

An example report in zip format can be found in the [supplement directory](https://github.com/kircherlab/cfDNA-UniFlow/tree/main/supplement). For viewing, the zip file needs to be extracted.

</div>
