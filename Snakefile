import os
import sys
from os import listdir
from os.path import isfile, join
import re
from Bio import SeqIO
import gzip
from collections import Counter

configfile: "config.yaml"
myfastqpath = "fastq/"
sys.path.append(config["ngs_path"]) # needed to correctly find helper script
import ngs_helper.ngs_helper as ngs_helper #provides the helper functions

# HELPER FUNCTIONS
# Create function for creating rule sets

def choose_rule_all(config):
    """
    Selects the input needed for 'rule all' in snakemake pipeline

    Input Parameter: 
    config (dict): Dictionary derived from config.yaml and any additional 
    key:value pairs added during the file preperation steps. 

    Returns:
    list: List of required input files for 'rule all'. 
    """
    myout = []
    if config["to_multiqc"] == "TRUE":
        myout.append("output/multiqc_report.html")
    if config["to_bw"] == "TRUE" and config["experiment"] != "rnaseq":
        myout.append(
            expand('output/bw/{sample}.unique.sorted.rmdup.chr.bw', 
                sample=SAMPLES))
    if config["to_bed"] == "TRUE" and config["experiment"] != "rnaseq":
        myout.append(
            expand('output/bed/{sample}.unique.sorted.rmdup.chr.bed', 
                sample=SAMPLES))
    if config["to_tdf"] == "TRUE" and config["experiment"] != "rnaseq":
        myout.append(
            expand('output/tdf/{sample}.unique.sorted.rmdup.tdf', 
                sample=SAMPLES))
    if config["experiment"] == "rnaseq":
        myout.append("output/counts_matrix.txt")
    return(myout)

# Create read-pair inputs for sample processing


def create_inputs(config):
    """
    Creates the fastq file inputs needed for read trimming steps of 
    the snakemake pipeline

    Input Parameter: 
    config (dict): Dictionary derived from config.yaml and any 
    additional key:value pairs added during the file preperation steps. 

    Returns:
    list: List of two strings; 
    1st string denotes the forward read
    2nd string denotes the reverse read
    """
    return([("fastq/{sample}" + expand("{ending}{suffix}", \
            ending=R1_file_ending, suffix=suffix)[0]+""),
            ("fastq/{sample}" + expand("{ending}{suffix}", \
            ending=R2_file_ending, suffix=suffix)[0]+"")])

# End helper functions


# Retrieve the list of fastq files
onlyfiles, gzfiles = ngs_helper.getfilelist(myfastqpath)

# Raise exception if no fastq files present
if len(gzfiles) == 0 and len(onlyfiles) == 0:
    raise NameError(
        "You do not seem to have any fastq files present to process. Exiting.")

# Raise exception if fastq files are a mixture of gzipped and non-gzipped files
if len(gzfiles) > 0 and len(gzfiles) != len(onlyfiles):
    myinput = "You have a mixture of gzipped files and non-gzipped files\n \
                Only {} of total {} files are gzipped!"
    raise NameError(print(myinput.format(len(gzfiles), len(onlyfiles))))

# Unify fastq file endings and return the final ending to be used.
if len(gzfiles) > 0:
    R1_file_ending, R2_file_ending, onlyfiles, gzfiles = \
            ngs_helper.fix_input_files(".gz", gzfiles, myfastqpath)
    suffix = ".gz"
else:
    R1_file_ending, R2_file_ending, onlyfiles, gzfiles = \
            ngs_helper.fix_input_files("", onlyfiles, myfastqpath)
    suffix = ""

sample_string = myfastqpath + "{sample}" + R1_file_ending + suffix
SAMPLES, = glob_wildcards(sample_string)

# Check the file pairing
# Raise exception for non-paired PE files
if config["type"] == "single":
    print("You selected single-end reads\nRead pairing not being checked...")
elif config["type"] == "paired":
    len_r1 = len([i for i in onlyfiles if i.endswith(R1_file_ending + suffix)])
    if len_r1*2 != len(onlyfiles):
        myinput = "One or more samples do not have a read pair!\nIf using \
            paired-end samples, please ensure each sample has read 1 and \
            read 2 files\nAborting..."
        raise NameError(myinput)  # Raise exception to break workflow
else:
    myinput = "You have specified unknown read type: " + \
        config["type"] + "\nPlease specify either \"paired\" or \"single\" \
        in the config.yaml file, then rerun the pipeline."
    raise NameError(myinput)

# Retrieve cut&run read lengths for use as parameters
if len(gzfiles) > 0 and config["experiment"] == "cutrun":
    config["read_length"], config["read_length_max"] = \
    ngs_helper.check_readlength(suffix, gzfiles, R1_file_ending, myfastqpath)
elif len(gzfiles) == 0 and config["experiment"] == "cutrun": 
    config["read_length"], config["read_length_max"] = \
    ngs_helper.check_readlength(suffix, onlyfiles, R1_file_ending, myfastqpath)

# Generate input rule for Snakemake
rule all:
    input:
        choose_rule_all(config)

config["suffix"] = suffix

if config["experiment"] == "cutrun":
    # Begin Snakemake pre-processing for Cut&Run samples
    if config["type"] == "paired":
        rule trim_fastq_fastqc:
            input:
                pair1 = (join(myfastqpath, "{sample}") + expand(
                    "{ending}{suffix}", ending=R1_file_ending, suffix=suffix)[0]),
                pair2 = (join(myfastqpath, "{sample}") + expand(
                    "{ending}{suffix}", ending=R2_file_ending, suffix=suffix)[0])
            output:
                trimmed_pair1 = temp(
                    "output/trim_fastq/{sample}_R1_val_1.fq.gz"),
                trimmed_pair2 = temp(
                    "output/trim_fastq/{sample}_R2_val_2.fq.gz"),
                fastqc_zipfile1 = "output/fastqc/{sample}_R1_fastqc.zip",
                fastqc_zipfile2 = "output/fastqc/{sample}_R2_fastqc.zip"
            log:
                "output/logs/{sample}.trim_adapters.log"
            params:
                config["suffix"]
            run:
                # mv files to R1 and R2 ending in temporary directory
                shell("mkdir -p output/temp_dir")
                shell("cp {input.pair1} \
                    output/temp_dir/{wildcards.sample}_R1.fq{params}")
                shell("cp {input.pair2} \
                    output/temp_dir/{wildcards.sample}_R2.fq{params}")
                shell("trim_galore \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{params} \
                    output/temp_dir/{wildcards.sample}_R2.fq{params} --paired \
                    -o ./output/trim_fastq")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{params} \
                    output/temp_dir/{wildcards.sample}_R2.fq{params} \
                    -o ./output/fastqc")

    if config["type"] == "single":
        rule trim_fastq_fastqc:
            input:
                pair1 = create_inputs(config)[0]
            output:
                trimmed_pair1 = temp(
                    "output/trim_fastq/{sample}_R1_val_1.fq.gz"),
                trimmed_pair2 = temp(
                    "output/trim_fastq/{sample}_R2_val_2.fq.gz"),
                fastqc_zipfile1 = "output/fastqc/{sample}_R1_fastqc.zip",
                fastqc_zipfile2 = "output/fastqc/{sample}_R2_fastqc.zip"
            log:
                "output/logs/{sample}.trim_adapters.log"
            params:
                config["suffix"]
            run:
                # mv files to R1 and R2 ending in temporary directory
                shell("mkdir -p output/temp_dir")
                shell("cp {input.pair1} \
                    output/temp_dir/{wildcards.sample}_R1.fq{params}")
                shell("trim_galore \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{params} \
                    -o ./output/trim_fastq --basename {wildcards.sample}")
                shell("mv output/trim_fastq/{wildcards.sample}_trimmed.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{params} \
                    -o ./output/fastqc")
                shell("touch {output.trimmed_pair2}")
                shell("touch {output.fastqc_zipfile2}")

    rule split_length_long:
        input:
            trimg_pair1 = "output/trim_fastq/{sample}_R1_val_1.fq.gz",
            trimg_pair2 = "output/trim_fastq/{sample}_R2_val_2.fq.gz"
        output:
            cut_r1_p1 = temp("output/trim_fastq/{sample}_t1_R1.len75.fastq"),
            cut_r2_p1 = temp("output/trim_fastq/{sample}_t1_R2.len75.fastq")
        params:
            read_length = config["read_length"],
            suffix = config["suffix"]
        log:
            "output/logs/{sample}.split_length_keeplong.log"
        run:
            #This step removes the temporary files in "output/temp_dir"
            if config["type"] == "paired":
                shell("cutadapt --minimum-length {params.read_length} \
                    -o {output.cut_r1_p1} {input.trimg_pair1}"),
                shell("cutadapt --minimum-length {params.read_length} \
                    -o {output.cut_r2_p1} {input.trimg_pair2}")
                shell(
                    "rm output/temp_dir/{wildcards.sample}_R1.fq{params.suffix}")
                shell(
                    "rm output/temp_dir/{wildcards.sample}_R2.fq{params.suffix}")
            else:
                shell("cutadapt --minimum-length {params.read_length} \
                    -o {output.cut_r1_p1} {input.trimg_pair1}")
                shell("touch {output.cut_r2_p1}")
                shell("rm output/fastqc/{wildcards.sample}_R2_fastqc.zip")
                shell(
                    "rm output/temp_dir/{wildcards.sample}_R1.fq{params.suffix}")

    rule split_length_short:
        input:
            trimg_pair1 = "output/trim_fastq/{sample}_R1_val_1.fq.gz",
            trimg_pair2 = "output/trim_fastq/{sample}_R2_val_2.fq.gz"
        output:
            cut_r1_p2 = temp("output/trim_fastq/{sample}_t1_R1.lt75.fastq"),
            cut_r2_p2 = temp("output/trim_fastq/{sample}_t1_R2.lt75.fastq")
        params:
            read_length = config["read_length_max"]
        log:
            "output/logs/{sample}.split_length_keepshort.log"
        run:
            if config["type"] == "paired":
                shell("cutadapt --maximum-length {params.read_length} \
                    -o {output.cut_r1_p2} {input.trimg_pair1}"),
                shell("cutadapt --maximum-length {params.read_length} \
                    -o {output.cut_r2_p2} {input.trimg_pair2}")
            else:
                shell("cutadapt --maximum-length {params.read_length} \
                    -o {output.cut_r1_p2} {input.trimg_pair1}")
                shell("touch {output.cut_r2_p2}")

    rule trim_long:
        input:
            cut_r1_p1 = "output/trim_fastq/{sample}_t1_R1.len75.fastq",
            cut_r2_p1 = "output/trim_fastq/{sample}_t1_R2.len75.fastq"
        output:
            cut_r1_p3 = temp(
                "output/trim_fastq/{sample}_t1_R1.len75_trim.fastq"),
            cut_r2_p3 = temp(
                "output/trim_fastq/{sample}_t1_R2.len75_trim.fastq")
        log:
            "output/logs/{sample}.trim_long.log"
        run:
            if config["type"] == "paired":
                shell("cutadapt -u -6 -o {output.cut_r1_p3} {input.cut_r1_p1}"),
                shell("cutadapt -u -6 -o {output.cut_r2_p3} {input.cut_r2_p1}")
            else:
                shell("cutadapt -u -6 -o {output.cut_r1_p3} {input.cut_r1_p1}"),
                shell("touch {output.cut_r2_p3}")

    rule combine_split_lengths:
        input:
            cut_r1_p3 = "output/trim_fastq/{sample}_t1_R1.len75_trim.fastq",
            cut_r2_p3 = "output/trim_fastq/{sample}_t1_R2.len75_trim.fastq",
            cut_r1_p2 = "output/trim_fastq/{sample}_t1_R1.lt75.fastq",
            cut_r2_p2 = "output/trim_fastq/{sample}_t1_R2.lt75.fastq"
        output:
            cut_r1_p4 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz",
            cut_r2_p4 = "output/trim_fastq/{sample}_R2_trimmed.fq.gz"
        run:
            if config["type"] == "paired":
                shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > \
                    output/trim_fastq/{wildcards.sample}_t2_R1.fastq"),
                shell("cat {input.cut_r2_p3} {input.cut_r2_p2} > \
                    output/trim_fastq/{wildcards.sample}_t2_R2.fastq"),
                shell("cat output/trim_fastq/{wildcards.sample}_t2_R1.fastq \
                    | paste - - - - > output/trim_fastq/{wildcards.sample}_t2_R1_flat.fastq"),
                shell("sort -k1,1 -T output/trim_fastq --parallel=20 -t \" \" \
                    output/trim_fastq/{wildcards.sample}_t2_R1_flat.fastq \
                    | tr \"\t\" \"\\n\" > output/trim_fastq/{wildcards.sample}_R1_trimmed.fq"),
                shell("rm output/trim_fastq/{wildcards.sample}_t2_R1_flat.fastq"),
                shell("cat output/trim_fastq/{wildcards.sample}_t2_R2.fastq \
                    | paste - - - - > output/trim_fastq/{wildcards.sample}_t2_R2_flat.fastq"),
                shell("sort -k1,1 -T output/trim_fastq --parallel=20 -t \" \" \
                    output/trim_fastq/{wildcards.sample}_t2_R2_flat.fastq \
                    | tr \"\t\" \"\\n\" > output/trim_fastq/{wildcards.sample}_R2_trimmed.fq"),
                shell("rm output/trim_fastq/{wildcards.sample}_t2_R2_flat.fastq"),
                shell(
                    "gzip output/trim_fastq/{wildcards.sample}_R1_trimmed.fq"),
                shell(
                    "gzip output/trim_fastq/{wildcards.sample}_R2_trimmed.fq"),
                shell("rm output/trim_fastq/{wildcards.sample}_t2_R1.fastq"),
                shell("rm output/trim_fastq/{wildcards.sample}_t2_R2.fastq")
            else:
                shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq")
                shell("gzip output/trim_fastq/{wildcards.sample}_R1_trimmed.fq")
                shell("touch {output.cut_r2_p4}")
    # end pre-processing for cut and run samples

    # Alignment for cut and run samples
    rule bowtie2:
        input:
            trimmed_pair1 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz",
            trimmed_pair2 = "output/trim_fastq/{sample}_R2_trimmed.fq.gz"
        params:
            index = config["bowtie2_index"]
        output:
            bam = temp("output/bam/{sample}.sorted.bam"),
            bambai = temp("output/bam/{sample}.sorted.bam.bai")
        threads:
            config["threads_for_alignment"]
        log:
            "output/logs/{sample}.alignment.log"
        run:
            if config["type"] == "paired":
                shell("bowtie2 -p {threads} --dovetail --phred33 \
                    -x {params.index} -1 {input.trimmed_pair1} \
                    -2 {input.trimmed_pair2} 2> {log} > \
                    output/bam/{wildcards.sample}.sam"),
                shell("samtools sort output/bam/{wildcards.sample}.sam \
                    | samtools view -bS - > output/bam/{wildcards.sample}.bam"),
                shell("rm output/bam/{wildcards.sample}.sam"),
                shell("samtools index output/bam/{wildcards.sample}.bam")
            else:
                shell("bowtie2 -p {threads} --dovetail --phred33 \
                    -x {params.index} -U {input.trimmed_pair1} \
                    2> {log} > output/bam/{wildcards.sample}.sam"),
                shell("samtools sort output/bam/{wildcards.sample}.sam \
                    | samtools view -bS - > output/bam/{wildcards.sample}.bam"),
                shell("rm output/bam/{wildcards.sample}.sam"),
                shell("samtools index output/bam/{wildcards.sample}.bam")
                shell("rm {input.trimmed_pair2}")
            if config["keep_fastq"] == "FALSE" and config["type"] == "paired":
                shell("rm {input.trimmed_pair1} {input.trimmed_pair2}")
            elif config["keep_fastq"] == "FALSE" and config["type"] == "single":
                shell("rm {input.trimmed_pair1}")
            # Now sorting the bam file
            shell("samtools view -bh -f 3 -F 4 -F 8 \
                output/bam/{wildcards.sample}.bam > \
                output/bam/{wildcards.sample}_mapped.bam")
            shell("samtools index output/bam/{wildcards.sample}_mapped.bam")
            shell("samtools sort output/bam/{wildcards.sample}_mapped.bam > \
                output/bam/{wildcards.sample}.sorted.bam")
            shell("samtools index output/bam/{wildcards.sample}.sorted.bam")
            shell("rm output/bam/{wildcards.sample}_mapped.bam*")
            if config["keep_unfiltered_bam"] == "FALSE":
                shell("rm output/bam/{wildcards.sample}.bam")
                shell("rm output/bam/{wildcards.sample}.bam.bai")

    # This rule removes duplicate alignments from Cut and Run samples
    rule rmdup:
        input:
            "output/bam/{sample}.sorted.bam"
        output:
            bam = "output/bam/{sample}.unique.sorted.rmdup.bam",
            bambai = "output/bam/{sample}.unique.sorted.rmdup.bam.bai"
        log:
            "output/logs/{sample}.rmdup.log"
        run:
            shell("samtools rmdup {input} {output.bam} 2> {log}")
            shell("samtools index {output.bam}")


# Pre-process RNA-seq or Chip-seq samples
if config["experiment"] == "rnaseq" or config["experiment"] == "chipseq":
    if config["type"] == "paired":
        # Trim adaptors from paired-end fastq files
        rule trim_fastq_fastqc:
            input:
                pair1 = create_inputs(config)[0],
                pair2 = create_inputs(config)[1]
            output:
                trimmed_pair1 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz",
                trimmed_pair2 = "output/trim_fastq/{sample}_R2_trimmed.fq.gz",
                fastqc_zipfile1 = "output/fastqc/{sample}_R1_fastqc.zip",
                fastqc_zipfile2 = "output/fastqc/{sample}_R2_fastqc.zip"
            log:
                "output/logs/{sample}.trim_adapters.log"
            params:
                config["suffix"]
            run:
                shell("mkdir -p output/temp_dir")
                shell("cp {input.pair1} \
                    output/temp_dir/{wildcards.sample}_R1.fq{params}")
                shell("cp {input.pair2} \
                    output/temp_dir/{wildcards.sample}_R2.fq{params}")
                shell("trim_galore \
                    output/temp_dir/{wildcards.sample}_R1.fq{params} \
                    output/temp_dir/{wildcards.sample}_R2.fq{params} --paired \
                    -o ./output/trim_fastq")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{params} \
                    output/temp_dir/{wildcards.sample}_R2.fq{params} -o ./output/fastqc")
                shell("mv output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz"),
                shell("mv output/trim_fastq/{wildcards.sample}_R2_val_2.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")

        # Align paired-end trimmed reads to genome
        rule fastq_to_bam:
            input:
                trimmed_pair1 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz",
                trimmed_pair2 = "output/trim_fastq/{sample}_R2_trimmed.fq.gz"
            params:
                index = config["index"],
                suffix = config["suffix"]
            output:
                bam = "output/bam/{sample}.bam",
                bambai = "output/bam/{sample}.bam.bai"
            threads: config["threads_for_alignment"]
            log:
                "output/logs/{sample}.alignment.log"
            run:
                shell("hisat2 -p {threads} -x {params.index} \
                    -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} \
                    -S output/bam/{wildcards.sample}.sam 2> {log}"),
                shell("samtools sort output/bam/{wildcards.sample}.sam | \
                    samtools view -bS - > {output.bam}"),
                shell("samtools index {output.bam}")
                shell("rm output/bam/{wildcards.sample}.sam")
                shell(
                    "rm output/temp_dir/{wildcards.sample}_R1.fq{params.suffix}")
                shell(
                    "rm output/temp_dir/{wildcards.sample}_R2.fq{params.suffix}")
                if config["keep_fastq"] == "FALSE":
                    shell("rm {input.trimmed_pair1} {input.trimmed_pair2}")


    elif config["type"] == "single":
        # Trim adaptors from single-end fastq files
        rule trim_fastq_fastqc:
            input:
                pair1 = create_inputs(config)[0]
            output:
                trimmed_pair1 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz",
                fastqc_zipfile1 = "output/fastqc/{sample}_R1_fastqc.zip"
            log:
                "output/logs/{sample}.trim_adapters.log"
            params:
                config["suffix"]
            run:
                shell("mkdir -p output/temp_dir")
                shell("cp {input.pair1} \
                    output/temp_dir/{wildcards.sample}_R1.fq{params}")
                shell("trim_galore \
                    output/temp_dir/{wildcards.sample}_R1.fq{params} \
                    -o ./output/trim_fastq"),
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{params} \
                    -o ./output/fastqc")

        # Align trimmed single-end reads to genome
        rule fastq_to_bam:
            input:
                trimmed_pair1 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz"
            params:
                index = config["index"]
            output:
                bam = "output/bam/{sample}.bam",
                bambai = "output/bam/{sample}.bam.bai"
            threads: config["threads_for_alignment"]
            log:
                "output/logs/{sample}.alignment.log"
            run:
                shell("hisat2 -p {threads} -x {params.index} \
                    -U {input.trimmed_pair1} \
                    -S output/bam/{wildcards.sample}.sam 2> {log}"),
                shell("samtools sort output/bam/{wildcards.sample}.sam | \
                    samtools view -bS - > {output.bam}"),
                shell("rm output/bam/{wildcards.sample}.sam"),
                shell("samtools index {output.bam}")
                shell(
                    "rm output/temp_dir/{wildcards.sample}_R1.fq{params.suffix}")
                if config["keep_fastq"] == "FALSE":
                    shell("rm {input.trimmed_pair1}")
    # Common pre-processing of ChiP-seq and RNA-seq reads complete

    if config["experiment"] == "chipseq":
        # Remove and sort the multimapped reads from ChIP-seq samples
        rule bam_to_unique_mapped:
            input:
                "output/bam/{sample}.bam"
            output:
                bam = temp("output/bam/{sample}.sorted.bam")
            run:
                if config["type"] == "paired":
                    # mapped with pair, unique
                    shell("samtools view -bh -f 3 -F 4 -F 8 -F 256 {input} > \
                        output/bam/{wildcards.sample}_filtered.bam")
                    shell("samtools sort -O BAM -o {output.bam} \
                        output/bam/{wildcards.sample}_filtered.bam")
                    shell("rm output/bam/{wildcards.sample}_filtered.bam")
                else:
                    # mapped with pair, unique locations
                    shell("samtools view -bh -F 4 -F 256 {input} > \
                        output/bam/{wildcards.sample}_filtered.bam")
                    shell("samtools sort -O BAM -o {output.bam} \
                        output/bam/{wildcards.sample}_filtered.bam")
                    shell("rm output/bam/{wildcards.sample}_filtered.bam")
                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm output/bam/{wildcards.sample}.bam \
                        output/bam/{wildcards.sample}.bam.bai")

        # Remove duplicate reads from ChIP-seq samples
        rule sortedbam_to_rmdup:
            input:
                sorted_bam = "output/bam/{sample}.sorted.bam"
            output:
                dup_removed = "output/bam/{sample}.unique.sorted.rmdup.bam"
            log:
                "output/logs/{sample}.rmdup.log"
            run:
                shell("samtools rmdup {input.sorted_bam} {output.dup_removed} \
                    2> {log}")
                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm -f {input.sorted_bam} {input.sorted_bam}.bai")

    elif config["experiment"] == "rnaseq" and config["type"] == "paired":
        # Remove duplicate reads from paired-end RNA-seq samples
        rule sortedbam_to_rmdup:
            input:
                "output/bam/{sample}.bam"
            output:
                "output/bam/{sample}.sorted.rmdup.bam"
            log:
                "output/logs/{sample}.rmdup.log"
            run:
                shell("samtools rmdup {input} {output} 2> {log}")
                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm -f {input} {input}.bai")

    elif config["experiment"] == "rnaseq" and config["type"] == "single":
        # Sort reads from single-end RNA-seq samples
        rule sortedbam_to_rmdup:
            input:
                "output/bam/{sample}.bam"
            output:
                "output/bam/{sample}.sorted.bam"
            log:
                "output/logs/{sample}.rmdup.log"
            run:
                shell("cp {input} {output}")
                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm -f {input} {input}.bai")

# Additional rules for ChIP-seq
# Create TDF files
rule rmdup_to_tdf:
    input:
        dup_removed = "output/bam/{sample}.unique.sorted.rmdup.bam"
    params:
        chr_sizes = config["chr_sizes"]
    output:
        tdf = "output/tdf/{sample}.unique.sorted.rmdup.tdf"
    log:
        "output/logs/{sample}.tdf.log"
    shell:
        "igvtools count {input.dup_removed} {output.tdf} {params.chr_sizes} "
        "2> {log}"

# Eliminate extraneous chromosomes from bam file
rule rmdup_to_chrbam:
    input:
        dup_removed = "output/bam/{sample}.unique.sorted.rmdup.bam"
    output:
        chrbam = "output/bam/{sample}.unique.sorted.rmdup.chr.bam",
        chrbambai = "output/bam/{sample}.unique.sorted.rmdup.chr.bam.bai"
    log:
        "output/logs/{sample}.chrbam.log"
    run:
        shell("samtools view -H {input.dup_removed} | \
            sed -e \"s/SN:\([0-9XY]\)/SN:chr\\1/\" -e \"s/SN:MT/SN:chrM/\" | \
            samtools reheader - {input.dup_removed} > {output.chrbam}")
        shell("samtools index {output.chrbam}")

# Create bigwig from bam file (main chromosomes only)
rule chrbam_to_bw:
    input:
        chrbam = "output/bam/{sample}.unique.sorted.rmdup.chr.bam",
        chrbambai = "output/bam/{sample}.unique.sorted.rmdup.chr.bam.bai"
    output:
        bw_file = "output/bw/{sample}.unique.sorted.rmdup.chr.bw"
    log:
        "output/logs/{sample}.bw.log"
    run:
        shell("bamCoverage -b {input.chrbam} -o {output.bw_file} --binSize 10 \
            --normalizeUsing RPKM")

# Create bed file from bam file (main chromosomes only)
rule chrbam_to_bed:
    input:
        chrbam = "output/bam/{sample}.unique.sorted.rmdup.chr.bam"
    output:
        bed = "output/bed/{sample}.unique.sorted.rmdup.chr.bed"
    log:
        "output/logs/{sample}.bed.log"
    shell:
        "bedtools bamtobed -i {input.chrbam} > {output.bed} 2> {log}"

# Create counts files for RNA-seq
config["gene_scheme"] = "-t gene"  # or "-g gene_id"
rule sortedbam_to_counts:
    input:
        sorted_bam = "output/bam/{sample}.sorted.bam" if config[
            "type"] == "single" else "output/bam/{sample}.sorted.rmdup.bam"
    output:
        counts = "output/counts/{sample}.counts.txt"
    params:
        gtf = config["gtf"],
        gene_scheme = config["gene_scheme"]
    log:
        "output/logs/{sample}.feature_counts.log"
    run:
        if config["count_scheme"] == "fraction" and config["type"] == "paired":
            shell("featureCounts -p -O --fraction {params.gene_scheme} \
                -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}")
        elif config["count_scheme"] == "fraction" and config["type"] == "single":
            shell("featureCounts -O --fraction {params.gene_scheme} \
                -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}")
        elif config["count_scheme"] == "count_all" and config["type"] == "paired":
            shell("featureCounts -p -O {params.gene_scheme} -a {params.gtf} \
                -o {output.counts} {input.sorted_bam} 2> {log}")
        elif config["count_scheme"] == "count_all" and config["type"] == "single":
            shell("featureCounts -O {params.gene_scheme} -a {params.gtf} \
                -o {output.counts} {input.sorted_bam} 2> {log}")
        elif config["count_scheme"] == "count_uniq" and config["type"] == "paired":
            shell("featureCounts -p {params.gene_scheme} -a {params.gtf} \
                -o {output.counts} {input.sorted_bam} 2> {log}")
        elif config["count_scheme"] == "count_uniq" and config["type"] == "single":
            shell("featureCounts {params.gene_scheme} -a {params.gtf} \
                -o {output.counts} {input.sorted_bam} 2> {log}")

# Compile counts for RNA-seq
rule counts_matrix:
    input:
        counts = expand("output/counts/{sample}.counts.txt", sample=SAMPLES)
    output:
        matrix = "output/counts_matrix.txt"
    params:
        config["gene_scheme"]
    run:
        import pandas as pd

        dict_of_counts = {}

        for file in input:
            sample = file.split(".")[0]
            dict_of_counts[sample] = {}

            with open(file, "r") as infile:
                next(infile)
                next(infile)
                for lines in infile:
                    lines = lines.strip().split("\t")
                    if {params} == "-t gene":
                        dict_of_counts[sample][lines[0]] = int(float(lines[7]))
                    else:
                        dict_of_counts[sample][lines[0]] = int(float(lines[6]))

        dataframe = pd.DataFrame(dict_of_counts)
        dataframe.to_csv(output[0], sep='\t')

# Create multiqc report for RNA-seq
if config["experiment"] == "rnaseq":
    rule run_multiqc:
        input:
            matrix = "output/counts_matrix.txt"
        output:
            multiqc_report = "output/multiqc_report.html"
        params:
            multiqc_config = config["multiqc_yaml"]
        shell:
            "multiqc . -f --config {params.multiqc_config}"

# Create multiqc report for ChIP-seq and Cut&Run samples
elif config["experiment"] == "chipseq" or config["experiment"] == "cutrun":
    rule run_multiqc:
        input:
            chrbam = expand(
                "output/bam/{sample}.unique.sorted.rmdup.chr.bam", sample=SAMPLES)
        output:
            multiqc_report = "output/multiqc_report.html"
        params:
            multiqc_config = config["multiqc_yaml"]
        shell:
            "multiqc . -f --config {params.multiqc_config}"
