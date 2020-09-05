import os
import sys
from os import listdir
from os.path import isfile, join, isdir
import re
from Bio import SeqIO
import gzip
from collections import Counter

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
        myout.append("output/counts_matrix.txt")
    if config["experiment"] != "rnaseq":
        myout.extend(expand("output/trim_fastq/{sample}_R{reads}_trimmed.fq.gz", sample=SAMPLES, reads=[1,2]))
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

## Determine if genome index should be built for STAR aligner
if config["use_star"] == "TRUE":
    #verify that star index exists
    if isfile(os.path.join(config["star_indexloc"],"genomeParameters.txt")):
        print("Using STAR index in "+config["star_indexloc"])
    else:
        # verify that fasta needed for genome does exist
        if isfile(config["genome_fasta"]):
            # star_index_needed="TRUE"        
            print("Generating STAR index from the following fasta file and GTF file:\n"+\
                config["genome_fasta"]+"\n"+config["gtf"])
            if isdir(config["star_indexloc"]):
                print("STAR index will be created in\n"+config["star_indexloc"])
            else:
                print("STAR index will be created in current directory.")
                config["star_indexloc"]="index/"
        else:
            print("Fasta file specified for creating STAR index was invalid. Now exiting...")
            raise NameError(print(config["genome_fasta"]))

##### RULES SECTION
# Generate input rule for Snakemake
rule all:
    input:
        choose_rule_all(config)

## Select correct rule(s) for trimming reads
if config["experiment"] == "cutrun":
    include: os.path.join(config["ngs_path"],"cr_rules.snake")
else:
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
            pair2 = create_fastq_inputs(config)[1]
        run:
            shell("mkdir -p output/temp_dir")
            if config["type"] == "paired":
                # mv files to R1 and R2 ending in temporary directory
                shell("cp {input.pair1} \
                    output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                shell("cp {params.pair2} \
                    output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
                shell("trim_galore \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    output/temp_dir/{wildcards.sample}_R2.fq{suffix} --paired \
                    -o ./output/trim_fastq")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    output/temp_dir/{wildcards.sample}_R2.fq{suffix} \
                    -o ./output/fastqc")
                shell("mv output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz"),
                shell("mv output/trim_fastq/{wildcards.sample}_R2_val_2.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")
            if config["type"] == "single":
                # mv files to R1 and R2 ending in temporary directory
                shell("cp {input.pair1} \
                    output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                shell("trim_galore \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/trim_fastq --basename {wildcards.sample}")
                shell("mv output/trim_fastq/{wildcards.sample}_trimmed.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/fastqc")
                shell("touch {output.trimmed_pair2}")
                shell("touch {output.fastqc_zipfile2}")

## Select correct rules for aligning reads
if config["use_star"] == "TRUE":
    # ruleorder: bam_to_unique_mapped > fastq_to_bam_STAR
    include: os.path.join(config["ngs_path"],"star_rules.snake")

if config["use_star"] == "FALSE":
    rule fastq_to_bam_HISAT:
        input:
            trimmed_pair = ["output/trim_fastq/{sample}_R1_trimmed.fq.gz", \
                            "output/trim_fastq/{sample}_R2_trimmed.fq.gz"]
        params:
            index = config["hisat2_index"]
        output:
            bam = "output/bam/{sample}.bam",
            bambai = "output/bam/{sample}.bam.bai"
        threads: config["threads_for_alignment"]
        log:
            "output/logs/{sample}.alignment.log"
        run:
            if config["experiment"] != "rnaseq" :
                # Splicing is not desired in cut&run and chipseq
                if config["type"] == "paired":
                    if config["cufflinks_bam"] == "FALSE" :
                        shell("hisat2 -p {threads} -x {params.index} \
                            -1 {input.trimmed_pair}[0] -2 {input.trimmed_pair}[1] \
                            --no-spliced-alignment \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                    else:
                        shell("hisat2 -p {threads} -x {params.index} \
                            --pen-noncansplice 1000000 \
                            --no-spliced-alignment \
                            -1 {input.trimmed_pair}[0] -2 {input.trimmed_pair}[1] \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                if config["type"] == "single":        
                    if config["cufflinks_bam"] == "FALSE":
                        shell("hisat2 -p {threads} -x {params.index} \
                            -U {input.trimmed_pair}[0] \
                            --no-spliced-alignment \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                    else:
                        shell("hisat2 -p {threads} -x {params.index} \
                            --pen-noncansplice 1000000 -U {input.trimmed_pair}[0] \
                            --no-spliced-alignment \
                            -S output/bam/{wildcards.sample}.sam 2> {log}") 
            if config["experiment"] == "rnaseq" :
                # Splicing is desired in rnaseq
                if config["type"] == "paired":
                    if config["cufflinks_bam"] == "FALSE" :
                        shell("hisat2 -p {threads} -x {params.index} \
                            -1 {input.trimmed_pair}[0] -2 {input.trimmed_pair}[1] \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                    else:
                        shell("hisat2 -p {threads} -x {params.index} \
                            --pen-noncansplice 1000000 \
                            -1 {input.trimmed_pair}[0] -2 {input.trimmed_pair}[1] \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                if config["type"] == "single":        
                    if config["cufflinks_bam"] == "FALSE":
                        shell("hisat2 -p {threads} -x {params.index} \
                            -U {input.trimmed_pair}[0] \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                    else:
                        shell("hisat2 -p {threads} -x {params.index} \
                            --pen-noncansplice 1000000 -U {input.trimmed_pair}[0] \
                            -S output/bam/{wildcards.sample}.sam 2> {log}") 
            shell("samtools sort output/bam/{wildcards.sample}.sam | \
                samtools view -bS - > {output.bam}"),
            shell("samtools index {output.bam}")
            shell("rm output/bam/{wildcards.sample}.sam")
            shell(
                "rm output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
            shell(
                "rm output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
            if config["keep_fastq"] == "FALSE":
                shell("rm {input.trimmed_pair}[0] {input.trimmed_pair}[1]")

# Remove and sort the multimapped reads from ChIP-seq and cut&run samples
if config["experiment"] != "rnaseq" :
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

    # Remove duplicate reads 
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

# Remove and sort the multimapped reads from RNA-seq samples
if config["experiment"] == "rnaseq":
    rule sortedbam_to_rmdup:
        input:
            "output/bam/{sample}.bam"
        output:
            "output/bam/{sample}.sorted.bam" if config["type"] == "single" else "output/bam/{sample}.sorted.rmdup.bam"
        log:
            "output/logs/{sample}.rmdup.log"
        run:
            if config["type"] == "paired":
                shell("samtools rmdup {input} {output} 2> {log}")
                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm -f {input} {input}.bai")
            else:
                shell("cp {input} {output}")
                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm -f {input} {input}.bai")            

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
        sorted_bam = "output/bam/{sample}.sorted.bam" if config[
            "type"] == "single" else "output/bam/{sample}.sorted.rmdup.bam"
    output:
        counts = "output/counts/{sample}.counts.txt"
    params:
        gtf = config["gtf"],
        gene_scheme = config["gene_scheme"] # for example, "gene" or "exon"
    log:
        "output/logs/{sample}.feature_counts.log"
    run:
        if config["count_scheme"] == "fraction":
            if config["type"] == "paired":
                shell("featureCounts -p -O --fraction  -t {params.gene_scheme} \
                    -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}")
            if config["type"] == "single":
                shell("featureCounts -O --fraction -t {params.gene_scheme} \
                    -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}")
        elif config["count_scheme"] == "all_reads":
            if config["type"] == "paired":
                shell("featureCounts -p -O -t {params.gene_scheme} -a {params.gtf} \
                    -o {output.counts} {input.sorted_bam} 2> {log}")
            if config["type"] == "single":
                shell("featureCounts -O -t {params.gene_scheme} -a {params.gtf} \
                    -o {output.counts} {input.sorted_bam} 2> {log}")
        elif config["count_scheme"] == "unique_reads":
            if config["type"] == "paired":
                shell("featureCounts -p -t {params.gene_scheme} -a {params.gtf} \
                    -o {output.counts} {input.sorted_bam} 2> {log}")
            if config["type"] == "single":
                shell("featureCounts -t {params.gene_scheme} -a {params.gtf} \
                    -o {output.counts} {input.sorted_bam} 2> {log}")

# Compile counts for RNA-seq
rule counts_matrix:
    input:
        counts = expand("output/counts/{sample}.counts.txt", sample=SAMPLES)
    output:
        matrix = "output/counts_matrix.txt"
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
                    dict_of_counts[sample][lines[0]] = int(float(lines[-1]))

        dataframe = pd.DataFrame(dict_of_counts)
        dataframe.to_csv(output[0], sep='\t')

# Create multiqc report (used for all workflows) 
rule run_multiqc:
    input:
        "output/counts_matrix.txt" if config["experiment"] == "rnaseq" else \
        expand("output/bam/{sample}.unique.sorted.rmdup.chr.bam", sample=SAMPLES)
    output:
        multiqc_report = "output/multiqc_report.html"
    params:
        multiqc_config = os.path.join(expand("{param}", param=config["ngs_path"])[0],"multiqc_config_template.yaml")
    shell:
        "multiqc . -f --outdir ./output/ --config {params.multiqc_config}"


