# RNAseq-Charmer:

This repository hosts an automated RNA-seq pipeline created using Snakemake. It can run on both cluster and a local machine.

## Dependency:
- [Anaconda](https://conda.io/docs/user-guide/install/linux.html) 

## Installation:
Clone this repository and change into the cloned RNAseq-Charmer directory. 

To create an environment using the environment.yaml file, type the following:

`conda env create -f environment.yaml`

This will create a conda environment called rnaseq_charmer.

## Usage on a local machine:

Copy the config.yaml, run_snakemake.sh and Snakefile to your RNA-seq project directory. This directory should also contain a directory called 'fastq' wherein all the fastq files have been placed. Make the required changes to the config.yaml file. 

Next, simply type `sh run_snakemake.sh` or `nohup sh run_snakemake.sh &` (to run in background).

## Usage on an LSF cluser:

Copy the config.yaml, run_snakemake_cluster.sh, cluster.json and Snakefile to your RNA-seq project directory. Again, this directory should also contain a directory called 'fastq' wherein all fastq files are placed. Make the required changes to the config.yaml and cluster.json file.

Next, simply type `nohup sh run_snakemake_cluster.sh &` (to run in background).
 
## Additional Snakemake options:

You can also customize the run_snakemake.sh and run_snakemake_cluster.sh according to your own needs. You might want to change the number of cores snakemake uses. Or you might want to do a dryrun. To access the additional options available in snakemake, type

`source activate rnaseq_charmer`

followed by 

`snakemake --help`
