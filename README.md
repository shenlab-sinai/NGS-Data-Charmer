# NGS-Data-Charmer:

This repository hosts an automated NGS data analysis pipeline (ChIP-seq, Iso-seq and RNA-seq) created using Snakemake. 

## Dependency:
- [Anaconda](https://conda.io/docs/user-guide/install/linux.html) 

## Installation:
Clone this repository and change into the cloned NGS-Data-Charmer directory. 

To create an environment using the environment.yaml file, type the following:

`conda env create -f environment.yaml`

This will create a conda environment called ngs_data_charmer.

## Usage note:

You must manually activate the conda environment prior to running the sh files. Type the following to activate the environment:

`conda activate ngs_data_charmer`

The reason for this requirement is a failure of the conda environment to successfully activate from within a shell script.

## Usage on a local machine:

Copy the config.yaml, run\_snakemake.sh and Snakefile to your NGS project directory. This directory should also contain a directory called 'fastq' wherein all the fastq files are placed. Make sure the project directory structure is as follows:
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

Copy the config.yaml, run\_snakemake\_cluster.sh, cluster.json and Snakefile to your NGS project directory. This directory should also contain a directory called 'fastq' wherein all fastq files are placed. Make sure the project directory structure is as follows:
```
.
├── cluster.json
├── config.yaml
├── fastq
│   ├── negD1-WC-40_S2_L003_R1_001.fastq.gz
│   └── negD1-WC-40_S2_L003_R2_001.fastq.gz
├── run_snakemake_cluster.sh
└── Snakefile (required for all analyses)
```
Make the required changes to the config.yaml and cluster.json file. 

Optionally, to utilize multiple cores on the cluster enviroment for computationally heavy tasks (such as alignment), you may change the number of cores utilized in the cluster.json file. Following is an example requesting 4 compute cores on a single node, asking 12GB of memory per core -

```
    "__default__" :
    {  
        "queue"     : "premium",
        "allocation": "acc_Nestlerlab",
        "n"         : 4,
        "resources" : "\"rusage[mem=12000] span[ptile=4]\"",
        "jobname"      : "{rule}.{wildcards}",
        "output"    : "logs/{rule}.{wildcards}.o",
        "error"     : "logs/{rule}.{wildcards}.e",
        "walltime"    : "02:00"
    }

```
You can make the above changes either to the '\__default\__' object alone or to any of the individual objects in the cluster.json file.

Finally, type `nohup sh run_snakemake_cluster.sh &` (to run in background).

## Steps in RNA-seq pipeline:

 ![ScreenShot](/dag/dag_rnaseq.png)

## Steps in ChIP-seq pipeline:

 ![ScreenShot](/dag/dag_chipseq.png)

## Steps in Cut&Run pipeline

 ![ScreenShot](/dag/dag_cutrun.png)

## Output directory structure:
```
.
├── cluster.json
├── config.yaml
├── fastq
│   ├── negD1-WC-40_S2_L003_R1_001.fastq.gz
│   └── negD1-WC-40_S2_L003_R2_001.fastq.gz
├── run_snakemake_cluster.sh
├── Snakefile
└── output
    ├── counts_matrix.txt (RNA-seq)
    ├── multiqc_report.html
    ├── trim_fastq
    	├── negD1-WC-40_S2_L003_R1_trimmed.fq.gz (optional)
    	├── negD1-WC-40_S2_L003_R1.fq.gz_trimming_report.txt
    	├── negD1-WC-40_S2_L003_R2_trimmed.fq.gz (optional)
    	└── negD1-WC-40_S2_L003_R2.fq.gz_trimming_report.txt
    ├── logs
    ├── fastqc
    	├── negD1-WC-40_S2_L003_R1_fastqc.html
    	├── negD1-WC-40_S2_L003_R1_fastqc.zip
    	├── negD1-WC-40_S2_L003_R2_fastqc.html
    	└── negD1-WC-40_S2_L003_R2_fastqc.zip
    ├── bam
    	├── negD1-WC-40_S2_L003.bam (optional)
    	└── negD1-WC-40_S2_L003.(unique.)sorted.rmdup.bam
    ├── bw
    	└── negD1-WC-40_S2_L003.unique.sorted.rmdup.bw
    ├── tdf
    	└── negD1-WC-40_S2_L003.unique.sorted.rmdup.tdf
    └── counts
    	└── negD1-WC-40_S2_L003.counts.txt (RNA-seq)
```

## File naming requirements

It is best practice to ensure that sequencing files are properly paired (if paired-end sequencing was performed) and are named with consistent fastq file endings (e.g. ".R1.fq.gz" and ".R2.fq.gz") prior to analyzing sequencing reads. However, the Snakemake pipeline is robust to mixtures of fastq file endings. It does this by detecting the most common forward read file ending (e.g. "\*.R1.fq.gz"), then renaming any files that do not conform to the most common fastq file ending.

Please note that input fastq file names that fail to conform to any of the expected fastq filename endings (for example, "\_R1\_001.fastq.gz",".R1.fq.gz", "fq.gz", "fastq.gz", "fastq", and "fq" are all examples of allowed fastq file endings) will be ignored by Snakemake (e.g. files named "Treatment.gz" or "Treatment.txt" or "Treatment" will be ignored). Mixtures of ".gz" and non-gzipped *fastq* files are NOT allowed, and will result in an error. If you do have such a mixture, please gzip the uncompressed fastq files with the "gzip" command OR unzip any compressed fastq files prior to running the pipeline. 

## Optional output files

You may want to retain the trimmed fastq files or unfiltered bam files. This is possible by modifying the "output" options in the configuration file to "TRUE" or "FALSE". 

## Cut&Run note
For very large Cut&Run sequencing runs, the walltime for the 'combine_split_lengths' may need to be increased. The 'combine_split_lengths' step does utilize parallel processing (4 threads) in order to speed up the read sorting. However, if you find that the pipeline fails at the 'combine_split_lengths', try increasing the walltime in the cluster.json file. 

## Test dataset

A small example dataset, composed of two downsampled paired end Cut&Run samples, with 250K paired-end reads in each sample, are publically available for use: 
https://drive.google.com/drive/folders/1stoAAxEDL4dtJ4Vjm6y0FnnCH--dM_MZ?usp=sharing

While originally derived from Cut&Run sequencing, the sample dataset can also be run using the "rnaseq" and "chipseq" options.

## Additional Snakemake options:

You can also customize the run\_snakemake.sh and run\_snakemake_cluster.sh scripts according to your own needs. You might wish to change the number of cores snakemake uses. Or you might want to do a dryrun. To explore additional options available in snakemake, type:

`conda activate ngs_data_charmer`

followed by 

`snakemake --help`
