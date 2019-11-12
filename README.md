# NGS-Data-Charmer:

This repository hosts an automated NGS data analysis pipeline (ChIP-seq, Iso-seq and RNA-seq) created using Snakemake. 

## Dependency:
- [Anaconda](https://conda.io/docs/user-guide/install/linux.html) 

## Installation:
Clone this repository and change into the cloned NGS-Data-Charmer directory. 

To create an environment using the environment.yaml file, type the following:

`conda env create -f environment.yaml`

This will create a conda environment called ngs_data_charmer_copy.

## Usage note:

You must manually activate the conda environment prior to running the sh files. Type the following to activate the environment:

`conda activate ngs_data_charmer_copy`

The reason for this requirement is a failure of the conda environment to successfully activate from within a shell script.

## Usage on a local machine:

Copy the config.yaml, run_snakemake.sh and Snakefile to your NGS project directory. This directory should also contain a directory called 'fastq' wherein all the fastq files are placed. Make sure the project directory structure is as follows:
```
.
├── config.yaml
├── fastq
│   ├── D1-WC_S2_L003_R1_001.fastq.gz
│   └── D1-WC_S2_L003_R2_001.fastq.gz
├── run_snakemake.sh
└── Snakefile
```
Make the required changes to the config.yaml file.

Finally, type `sh run_snakemake.sh` followed by the maximum number of CPU cores to be used by snakemake. For example, type `sh run_snakemake.sh 2` for 2 CPU cores. You can also type `nohup sh run_snakemake.sh 2 &` to run the pipeline in background.

## Usage on an LSF cluser:

Copy the config.yaml, run_snakemake_cluster.sh, cluster.json and Snakefile to your NGS project directory. This directory should also contain a directory called 'fastq' wherein all fastq files are placed. Make sure the project directory structure is as follows:
```
.
├── cluster.json
├── config.yaml
├── fastq
│   ├── negD1-WC-40_S2_L003_R1_001.fastq.gz
│   └── negD1-WC-40_S2_L003_R2_001.fastq.gz
├── run_snakemake_cluster.sh
└── Snakefile
```
Make the required changes to the config.yaml and cluster.json file.

Next, type `nohup sh run_snakemake_cluster.sh &` (to run in background).

## Steps in RNA-seq pipeline:

 ![ScreenShot](/dag/dag_rnaseq.png)

## Steps in ChIP-seq pipeline:

 ![ScreenShot](/dag/dag_chipseq.png)

## Steps in Cut&Run pipeline

 ![ScreenShot](/dag/dag_cutrun.png)

## File naming requirements

 The Snakemake pipeline is now able to handle mixtures of fastq file endings. It does this by detecting the most common forward read file ending (e.g. *.R1.fq.gz), then renaming any files that do not conform to that file ending. This process allows for a large variety of possible fastq file endings, however, if you use a completely unexpected file naming system (e.g. Sample1.A.fq.gz and Sample1.B.fq.gz, where you intended "A" to mean the forward read and "B" to mean the reverse read), the unexpected file names will be treated as single-end files.  Please note that files ending only in ".fq.gz", ".fastq.gz", ".fq", and ".fastq" (e.g. "SampleA.fq.gz", "SampleB.fastq.gz", "SampleC.fq", and "SampleD.fastq") will be assumed to be single end and be renamed to the most common forward read file ending (or if your dataset consists solely of such files, they will be renamed to "*.fq" or "*.fq.gz")

Mixtures of ".gz" and non-gzipped fastq files are not allowed, and will result in an error. If you do have such a mixture, please gzip the uncompressed fastq files with the "gzip" command OR unzip any compressed fastq files prior to running the pipeline. 

Please note that input fastq file names that do not conform to any of the expected patterns will be ignored by Snakemake (e.g. "Treatment.gz" or "Treatment.txt" or "Treatment"). Mixtures 

## Output options

You may want to retain the trimmed fastq files or unfiltered bam files. Additionally, you may want to only run the workflow until a certain point (e.g. until duplicate removal). This is possible by modifying the "output" options in the configuration file to "TRUE" or "FALSE". 

## Additional Snakemake options:

You can also customize the run_snakemake.sh and run_snakemake_cluster.sh scripts according to your own needs. You might wish to change the number of cores snakemake uses. Or you might want to do a dryrun. To explore additional options available in snakemake, type:

`source activate ngs_data_charmer_copy`

followed by 

`snakemake --help`
