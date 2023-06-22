import os
import sys
from os import listdir
from os.path import isfile, join, isdir
import re
from Bio import SeqIO
import gzip
from collections import Counter
from datetime import datetime

configfile: "config.yaml"
myfastqpath = "fastq/"
sys.path.append(config["ngs_path"]) # needed to correctly find helper script
import ngs_helper.ngs_helper as ngs_helper #provides the helper functions

## HELPER FUNCTIONS
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
        # myout.append(expand("output/bam/{sample}.unique.sorted.rmdup.chr.bam", sample=SAMPLES))
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
        myout.append("output/counts_genic_matrix.txt")
        myout.append("output/counts_exonic_matrix.txt")
    # if config["experiment"] != "rnaseq":
    #     myout.extend(expand("output/trim_fastq/{sample}_R{reads}_trimmed.fq.gz", sample=SAMPLES, reads=[1,2]))
    return(myout)

# Create read-pair inputs for sample processing
def create_fastq_inputs(config):
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
    return([os.path.join("fastq","{sample}"+expand("{ending}{suffix}", \
        ending=R1_file_ending, suffix=suffix)[0]+""),
        os.path.join("fastq","{sample}"+expand("{ending}{suffix}", \
            ending=R2_file_ending, suffix=suffix)[0]+"")])
## End helper functions

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

sample_string = os.path.join("fastq","{sample}"+R1_file_ending+suffix)
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
    read_length, read_length_max = \
    ngs_helper.check_readlength(suffix, gzfiles, R1_file_ending, myfastqpath)
elif len(gzfiles) == 0 and config["experiment"] == "cutrun": 
    read_length, read_length_max = \
    ngs_helper.check_readlength(suffix, onlyfiles, R1_file_ending, myfastqpath)

## Set up the alignment sensitivity parameter
if config["use_very_sensitive"] == "TRUE":
    sensitivity_level = "--very-sensitive"
else:
    sensitivity_level = ""

##### RULES SECTION
# Generate input rule for Snakemake
rule all:
    input:
        choose_rule_all(config)

## Select correct rule(s) for trimming reads
rule trim_fastq_fastqc:
    input:
        pair1 = create_fastq_inputs(config)[0]
    output:
        trimmed_pair1 = temp("output/trim_fastq/{sample}_R1_trimmed.fq.gz"),
        trimmed_pair2 = temp("output/trim_fastq/{sample}_R2_trimmed.fq.gz"),
        fastqc_zipfile1 = "output/fastqc/{sample}_R1_fastqc.zip",
        fastqc_zipfile2 = "output/fastqc/{sample}_R2_fastqc.zip"
    log:
        "output/logs/{sample}.trim_adapters.log"
    params:
        pair2 = create_fastq_inputs(config)[1],
        umi_1 = config["UMI_read1_pattern"],
        umi_2 = config["UMI_read2_pattern"]
    run:
        shell("mkdir -p output/temp_dir")
        if config['trim_polyA'] == "TRUE":
            if config["type"] == "paired":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools extract -I {input.pair1} --bc-pattern={params.umi_1}  --bc-pattern2={params.umi_2} \
                    --read2-in={params.pair2} --stdout=output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    --read2-out=output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
                else:
                    # mv files to R1 and R2 ending in temporary directory
                    shell("cp {input.pair1} \
                        output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                    shell("cp {params.pair2} \
                        output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
                shell("trim_galore \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    output/temp_dir/{wildcards.sample}_R2.fq{suffix} --paired \
                    -o ./output/trim_fastq")
                shell("trim_galore \
                    ./output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz \
                    ./output/trim_fastq/{wildcards.sample}_R2_val_2.fq.gz --paired --polyA \
                     --basename {wildcards.sample}_pat")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    output/temp_dir/{wildcards.sample}_R2.fq{suffix} \
                    -o ./output/fastqc")
                shell("mv ./{wildcards.sample}_pat_val_1.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz"),
                shell("mv ./{wildcards.sample}_pat_val_2.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")
                shell("mv {wildcards.sample}_R2_val_2.fq.gz_trimming_report.txt ./output/trim_fastq/")
                shell("mv {wildcards.sample}_R1_val_1.fq.gz_trimming_report.txt ./output/trim_fastq/")
                shell("rm ./output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz")
                shell("rm ./output/trim_fastq/{wildcards.sample}_R2_val_2.fq.gz")
            if config["type"] == "single":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools extract --stdin={input.pair1} --bc-pattern={params.umi_1} \
                        --log={log} --stdout output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                else:
                    # mv files to R1 and R2 ending in temporary directory
                    shell("cp {input.pair1} \
                        output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                shell("trim_galore \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/trim_fastq --basename {wildcards.sample}")
                shell("trim_galore --polyA \
                    ./output/trim_fastq/{wildcards.sample}_trimmed.fq.gz ") # new rule --basename {wildcards.sample}_pat
                shell("mv {wildcards.sample}_trimmed_trimmed.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz")
                shell("mv {wildcards.sample}_trimmed.fq.gz_trimming_report.txt ./output/trim_fastq/")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/fastqc")
                shell("touch {output.trimmed_pair2}")
                shell("touch {output.fastqc_zipfile2}")
        else:
            if config["type"] == "paired":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools extract -I {input.pair1} --bc-pattern={params.umi_1}  --bc-pattern2={params.umi_2} \
                    --read2-in={params.pair2} --stdout=output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    --read2-out=output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
                else:
                    # mv files to R1 and R2 ending in temporary directory
                    shell("cp {input.pair1} \
                        output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                    shell("cp {params.pair2} \
                        output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
                if config["experiment"] == "cutrun":
                    shell("trim_galore --clip_R1 6 --clip_R2 6 \
                        --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                        output/temp_dir/{wildcards.sample}_R2.fq{suffix} --paired --trim-n \
                        -o ./output/trim_fastq")
                else:
                    shell("trim_galore \
                        --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                        output/temp_dir/{wildcards.sample}_R2.fq{suffix} --paired --trim-n \
                        -o ./output/trim_fastq")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    output/temp_dir/{wildcards.sample}_R2.fq{suffix} \
                    -o ./output/fastqc")
                shell("mv output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz"),
                shell("mv output/trim_fastq/{wildcards.sample}_R2_val_2.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")
            if config["type"] == "single":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools extract --stdin={input.pair1} --bc-pattern={params.umi_1} \
                        --log={log} --stdout output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                else:
                    # mv files to R1 and R2 ending in temporary directory
                    shell("cp {input.pair1} \
                        output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                shell("trim_galore --trim-n \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/trim_fastq --basename {wildcards.sample}")
                shell("mv output/trim_fastq/{wildcards.sample}_trimmed.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/fastqc")
                shell("touch {output.trimmed_pair2}")
                shell("touch {output.fastqc_zipfile2}")


rule fastq_to_bam_HISAT:
    input:
        trimmed_pair = ["output/trim_fastq/{sample}_R1_trimmed.fq.gz", \
                            "output/trim_fastq/{sample}_R2_trimmed.fq.gz"]
    params:
        hisat2index = config["hisat2_index"]
    output:
        bam = "output/bam/{sample}.bam",
        bambai = "output/bam/{sample}.bam.bai"
    threads: config["threads_for_alignment"]
    log:
        "output/logs/{sample}.alignment.log"
    run:
        timelog="output/logs/"+wildcards.sample+".alignment.log"
        # print(timelog)
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        starttime = "Alignment starting date and time =" + dt_string
        ## Timestamp in log file, using datetime object containing current date and time
        ## ChIP-Seq alignment rules
        if config["experiment"] != "rnaseq" :
            if config["type"] == "paired":
                # Perform the alignment
                # Splicing is not desired in cut&run and chipseq
                if config["cufflinks_bam"] == "FALSE" :
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            -1 {input.trimmed_pair[0]} \
                            -2 {input.trimmed_pair[1]} \
                            --no-spliced-alignment \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                else:
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            --pen-noncansplice 1000000 \
                            --no-spliced-alignment \
                            -1 {input.trimmed_pair[0]} \
                            -2 {input.trimmed_pair[1]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
            if config["type"] == "single":        
                if config["cufflinks_bam"] == "FALSE":
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            -U {input.trimmed_pair[0]} \
                            --no-spliced-alignment \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                else:
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            --pen-noncansplice 1000000 -U {input.trimmed_pair[0]} \
                            --no-spliced-alignment \
                            -S output/bam/{wildcards.sample}.sam 2> {log}") 
        ## RNA-seq alignment rules
        if config["experiment"] == "rnaseq" :
            # Splicing is desired in rnaseq
            if config["type"] == "paired":
                if config["cufflinks_bam"] == "FALSE" :
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            -1 {input.trimmed_pair[0]} -2 {input.trimmed_pair[1]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                else:
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            --pen-noncansplice 1000000 \
                            -1 {input.trimmed_pair[0]} -2 {input.trimmed_pair[1]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
            if config["type"] == "single":        
                if config["cufflinks_bam"] == "FALSE":
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            -U {input.trimmed_pair[0]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                else:
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            --pen-noncansplice 1000000 -U {input.trimmed_pair[0]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
        # Timestamp in log file
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        endtime = "Alignment ending date and time =" + dt_string
        with open(timelog, 'a') as f:
            print(starttime, file=f)
            print(endtime, file=f)
        #
        ## Convert and cleanup the alignment files                             
        shell("samtools sort -@ 8 -O BAM -o {output.bam} output/bam/{wildcards.sample}.sam")
        shell("rm output/bam/{wildcards.sample}.sam")
        shell("samtools index {output.bam}")
        if config["experiment"] != "cutrun" :
            shell("rm output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
            if config["type"] == "paired":
                shell("rm output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
        if config["keep_fastq"] == "FALSE":
            shell("rm output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")

# Remove and sort the multimapped reads from ChIP-seq and cut&run samples
if config["experiment"] != "rnaseq" :
    rule bam_to_unique_mapped:
        input:
            "output/bam/{sample}.bam"
        output:
            bam = temp("output/bam/{sample}.sorted.bam")
        run:
            if config["type"] == "paired":
                if config["experiment"] == "chipseq" :
	            # mapped with pair, unique
                    shell("samtools view -bh -f 3 -F 4 -F 8 -F 256 {input} > \
                       output/bam/{wildcards.sample}_filtered.bam")
                    shell("samtools sort -O BAM -o {output.bam} \
                    	output/bam/{wildcards.sample}_filtered.bam")
	            shell("rm output/bam/{wildcards.sample}_filtered.bam")
                if config["experiment"] == "cutrun" :
                    ## Split file mapped pairs into discordant/cordant components
                    shell("samtools view -bh -f 3 -F 4 -F 8 -F 256 {input} > \
                        output/bam/{wildcards.sample}_mapped_paired.bam")
                    shell("samtools view -bh -f 1 -F 2 -F 4 -F 8 -F 256 {input} > \
                        output/bam/{wildcards.sample}_mapped_dis.bam")
                    ## remove overlapping bases from discordant pairs
                    shell("bam clipOverlap --in output/bam/{wildcards.sample}_mapped_dis.bam --out \
                        output/bam/{wildcards.sample}_clipped_dis.bam --stats --overlapsOnly")
                    ## combine non-discord and filtered discordant pairs
                    shell("samtools merge -f -@ 4 output/bam/{wildcards.sample}_filtered.bam \
                                output/bam/{wildcards.sample}_mapped_paired.bam \
                                output/bam/{wildcards.sample}_clipped_dis.bam")
                    ## sort combined bam file
                    shell("samtools sort -O BAM -o {output.bam} \
                        output/bam/{wildcards.sample}_filtered.bam")
                    shell("rm output/bam/{wildcards.sample}_filtered.bam \
                        output/bam/{wildcards.sample}_mapped_paired.bam \
                        output/bam/{wildcards.sample}_mapped_dis.bam \
                        output/bam/{wildcards.sample}_clipped_dis.bam")
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

    # Remove duplicate reads 
    rule sortedbam_to_rmdup:
        input:
            sorted_bam = "output/bam/{sample}.sorted.bam"
        output:
            dup_removed = "output/bam/{sample}.unique.sorted.rmdup.bam"
        log:
            "output/logs/{sample}.rmdup.log"
        run:
            if config["experiment"] == "cutrun" :
                shell("picard MarkDuplicates --REMOVE_DUPLICATES true -I {input.sorted_bam} -O {output.dup_removed} \
                    -M output/logs/{wildcards.sample}_marked_dup_metrics.txt -VALIDATION_STRINGENCY SILENT 2> {log}")
            if config["experiment"] == "chipseq" :
                shell("picard MarkDuplicates --REMOVE_DUPLICATES true -I {input.sorted_bam} -O {output.dup_removed} \
                    -M output/logs/{wildcards.sample}_marked_dup_metrics.txt 2> {log}")
            if config["keep_unfiltered_bam"] == "FALSE":
                shell("rm -f {input.sorted_bam} {input.sorted_bam}.bai")
            ## Rule for calculating insert size for PE sequencing
            if config["type"] == "paired":
                shell("picard CollectInsertSizeMetrics \
                    I=output/bam/{wildcards.sample}.unique.sorted.rmdup.bam \
                    O=output/logs/{wildcards.sample}_insert_size_metrics.txt \
                    H=output/logs/{wildcards.sample}_insert_size_histogram.pdf \
                    M=0.05 VALIDATION_STRINGENCY=SILENT")


# Remove and sort the multimapped reads from RNA-seq samples
if config["experiment"] == "rnaseq":
    rule sortedbam_to_rmdup:
        input:
            "output/bam/{sample}.bam"
        output:
            "output/bam/{sample}.sorted.bam" if config["type"] == "single" else "output/bam/{sample}.sorted.rmdup.bam"
        log:
            "output/logs/{sample}.rmdup.log"
        params:
            refflat = config["refflat"]
        run:
            if config["type"] == "paired":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools dedup --stdin={input} --log={log} --paired --unmapped-reads=use > {output}")
                else:
                    # shell("samtools markdup -r {input} {output} 2> {log}")
                    shell("picard MarkDuplicates REMOVE_DUPLICATES=true I={input} O={output} \
                        M=output/logs/{wildcards.sample}_marked_dup_metrics.txt 2> {log}")

                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm -f {input} {input}.bai")
                ## Rule for calculating insert size for PE sequencing
                shell("picard CollectInsertSizeMetrics \
                    I=output/bam/{wildcards.sample}.sorted.rmdup.bam \
                    O=output/logs/{wildcards.sample}_insert_size_metrics.txt \
                    H=output/logs/{wildcards.sample}_insert_size_histogram.pdf \
                    M=0.05")
                if config["refflat"] != "FALSE":
                    shell("picard CollectRnaSeqMetrics \
                    I=output/bam/{wildcards.sample}.sorted.rmdup.bam \
                    O=output/logs/{wildcards.sample}.RNA_Metrics \
                    REF_FLAT={params.refflat} \
                    STRAND=SECOND_READ_TRANSCRIPTION_STRAND")
            else:
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools dedup --stdin={input} --log={log} > {output}")
                else:
                    shell("cp {input} {output}")
                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm -f {input} {input}.bai")
                if config["refflat"] != "FALSE":
                    shell("picard CollectRnaSeqMetrics \
                    I=output/bam/{wildcards.sample}.sorted.bam \
                    O=output/logs/{wildcards.sample}.RNA_Metrics \
                    REF_FLAT={params.refflat} \
                    STRAND=FIRST_READ_TRANSCRIPTION_STRAND")

# Additional rules for ChIP-seq
# Create TDF files
rule rmdup_to_tdf:
    input:
        "output/bam/{sample}.unique.sorted.rmdup.bam"
    params:
        chr_sizes = config["chr_sizes"]
    output:
        tdf = "output/tdf/{sample}.unique.sorted.rmdup.tdf"
    log:
        "output/logs/{sample}.tdf.log"
    shell:
        "igvtools count {input} {output.tdf} {params.chr_sizes} "
        "2> {log}"

# Eliminate extraneous chromosomes from bam file
rule rmdup_to_chrbam:
    input:
        "output/bam/{sample}.unique.sorted.rmdup.bam"
    output:
        chrbam = "output/bam/{sample}.unique.sorted.rmdup.chr.bam",
        chrbambai = "output/bam/{sample}.unique.sorted.rmdup.chr.bam.bai"
    log:
        "output/logs/{sample}.chrbam.log"
    run:
        shell("samtools view -H {input} | \
            sed -e \"s/SN:\([0-9XY]\)/SN:chr\\1/\" -e \"s/SN:MT/SN:chrM/\" | \
            samtools reheader - {input} > {output.chrbam}")
        shell("samtools index {output.chrbam}")

# Create bigwig from bam file (main chromosomes only)
rule chrbam_to_bw:
    input:
        chrbam = rules.rmdup_to_chrbam.output.chrbam,
        chrbambai = rules.rmdup_to_chrbam.output.chrbambai
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

## Additional rules for RNA-seq
# Create counts files for RNA-seq
rule sortedbam_to_counts:
    input:
        sorted_bam = "output/bam/{sample}.sorted.bam" if config["type"] == "single" else "output/bam/{sample}.sorted.rmdup.bam"
    output:
        gene_counts = "output/counts/{sample}.gene.counts.txt",
        exon_counts = "output/counts/{sample}.exon.counts.txt"
    params:
        gtf = config["gtf"]
    log:
        "output/logs/{sample}.feature_counts.log"
    run:
        if config["count_scheme"] == "fraction":
            if config["type"] == "paired":
                shell("featureCounts -p -O --fraction  -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -p -O --fraction  -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
            if config["type"] == "single":
                shell("featureCounts -O --fraction -t gene \
                    -a {params.gtf} -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -O --fraction -t exon \
                    -a {params.gtf} -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
        elif config["count_scheme"] == "all_reads":
            if config["type"] == "paired":
                shell("featureCounts -p -O -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -p -O -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
            if config["type"] == "single":
                shell("featureCounts -O -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -O -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
        elif config["count_scheme"] == "unique_reads":
            if config["type"] == "paired":
                shell("featureCounts -p -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -p -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
            if config["type"] == "single":
                shell("featureCounts -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")

# Compile counts for RNA-seq
rule counts_matrix:
    input:
        gene_counts = expand("output/counts/{sample}.gene.counts.txt", sample=SAMPLES),
        exon_counts = expand("output/counts/{sample}.exon.counts.txt", sample=SAMPLES)
    output:
        gene_matrix = "output/counts_genic_matrix.txt",
        exon_matrix = "output/counts_exonic_matrix.txt"
    run:
        import pandas as pd
        import platform

        ## Process the genic counts
        dict_of_counts = {}
        for file in input.gene_counts:
            sample = file.split(".")[0]
            if platform.system() != 'Windows':
                sample = sample.split("/")[2]
            else:
                sample = sample.split("\\")[2]
            
            dict_of_counts[sample] = {}
            with open(file, "r") as infile:
                next(infile)
                next(infile)
                for lines in infile:
                    lines = lines.strip().split("\t")
                    dict_of_counts[sample][lines[0]] = int(float(lines[-1]))

        dataframe = pd.DataFrame(dict_of_counts)
        dataframe.to_csv(output[0], sep='\t')

        ## Process the exonic counts
        dict_of_counts = {}
        for file in input.exon_counts:
            sample = file.split(".")[0]
            if platform.system() != 'Windows':
                sample = sample.split("/")[2]
            else:
                sample = sample.split("\\")[2]
            
            dict_of_counts[sample] = {}
            with open(file, "r") as infile:
                next(infile)
                next(infile)
                for lines in infile:
                    lines = lines.strip().split("\t")
                    dict_of_counts[sample][lines[0]] = int(float(lines[-1]))

        dataframe = pd.DataFrame(dict_of_counts)
        dataframe.to_csv(output[1], sep='\t')

rule fpkm_matrix:
    input:
        gene_counts = "output/counts_genic_matrix.txt",
        exon_counts = "output/counts_exonic_matrix.txt",
        length_gene = expand("output/counts/{sample}.gene.counts.txt", sample=SAMPLES),
        length_exon = expand("output/counts/{sample}.gene.counts.txt", sample=SAMPLES)
    output:
        gene_fpkm = "output/fpkm_genic_matrix.txt",
        exon_fpkm = "output/fpkm_exonic_matrix.txt"
    run:
        import pandas as pd
        counts = pd.read_csv(input.gene_counts, sep='\t')
        gl = pd.read_csv(input.length_gene[0], comment='#', header=0, sep='\t')
        counts.sort_values('Unnamed: 0',axis=0,inplace=True)
        gl.sort_values('Geneid',axis=0,inplace=True)

        if gl.iloc[:,0].tolist()==counts.iloc[:,0].tolist():
            genelength = gl["Length"].values/1000
            counts_rev = counts.drop("Unnamed: 0", axis=1).copy() # remove identifier column
            colsum = counts_rev.copy().sum()/(10**6) # Sum count columns
            rpkm = counts_rev.divide(genelength,axis=0).divide(colsum,axis=1)
            dataframe = pd.concat([gl.iloc[:,0],rpkm],axis=1)
            dataframe.to_csv(output.gene_fpkm, index=False, sep='\t')
        else:
            print("Gene IDs not identically aligned\nFPKM cannot be generated.")

        counts = pd.read_csv(input.exon_counts, sep='\t')
        gl = pd.read_csv(input.length_exon[0], comment='#', header=0, sep='\t')
        counts.sort_values('Unnamed: 0',axis=0,inplace=True)
        gl.sort_values('Geneid',axis=0,inplace=True)

        if gl.iloc[:,0].tolist()==counts.iloc[:,0].tolist():
            genelength = gl["Length"].values/1000
            counts_rev = counts.drop("Unnamed: 0", axis=1).copy() # remove identifier column
            colsum = counts_rev.copy().sum()/(10**6) # Sum count columns
            rpkm = counts_rev.divide(genelength,axis=0).divide(colsum,axis=1)
            dataframe = pd.concat([gl.iloc[:,0],rpkm],axis=1)
            dataframe.to_csv(output.exon_fpkm, index=False, sep='\t')
        else:
            print("Gene IDs not identically aligned\nFPKM cannot be generated.")

# Create multiqc report (used for all workflows) 
rule run_multiqc:
    input:
        "output/fpkm_genic_matrix.txt" if config["experiment"] == "rnaseq" else \
        expand("output/bam/{sample}.unique.sorted.rmdup.chr.bam", sample=SAMPLES)
    output:
        multiqc_report = "output/multiqc_report.html"
    params:
        multiqc_config = os.path.join(expand("{param}", param=config["ngs_path"])[0],"multiqc_config_template.yaml")
    shell:
        "multiqc . -f --outdir ./ --config {params.multiqc_config}"


