# RNAseq-Charmer:
==============

This repository hosts an automated RNA-seq pipeline created using Snakemake. This pipeline currently runs on a local machine.

## Dependency:
- [Anaconda](https://conda.io/docs/user-guide/install/linux.html) 

## Installation:
To start using the pipeline, click on the green Clone or download icon at the top right corner of this page, and copy the web URL. On your terminal, navigate to the directory where you would like to place this pipeline and type 'git clone' followed by the copied URL. Change into the RNAseq-Charmer directory. 

To create an environment from the environment.yaml file, type the following:

'conda env create -f environment.yaml'

This will create a conda environment called rnaseq_charmer. Next, place the config.yaml and Snakefile into your project directory. The project directory should also contain the fastq files within a directory called 'fastq'.





