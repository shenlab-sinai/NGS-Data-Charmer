# RNAseq-Charmer:

This repository hosts an automated RNA-seq pipeline created using Snakemake. The pipeline currently runs on a local machine, and supports single end data only.

## Dependency:
- [Anaconda](https://conda.io/docs/user-guide/install/linux.html) 

## Installation:
Clone this repository onto your local desktop and change into the cloned RNAseq-Charmer directory. 

To create an environment using the environment.yaml file, type the following:

`conda env create -f environment.yaml`

This will create a conda environment called rnaseq_charmer.

## Usage:

Copy the config.yaml and Snakefile to your RNA-seq project directory. This directory should also contain a directory called 'fastq' wherein all the fastq files have been placed. Make the required changes to the config.yaml file. Next, to activate the environment 'rnaseq_charmer', type the following:

`source activate rnaseq_charmer`

Once the environment is activated, you are now ready to run the pipeline! Simply type `snakemake` or `nohup snakemake &` (to run in background). 
