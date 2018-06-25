SAMPLES, = glob_wildcards("fastq/{sample}_R1_001.fastq.gz")
configfile: "config.yaml"

ALL_TDF = expand('processed/{sample}.unique.sorted.rmdup.tdf', sample=SAMPLES)
ALL_BED = expand('processed/{sample}.unique.sorted.rmdup.chr.bed', sample=SAMPLES)

if config["experiment"] == "chipseq": # defining target
	rule all:
		input: ALL_BED + ALL_TDF

elif config["experiment"] == "rnaseq":
	rule all:
		input:
			"processed/htseq_counts_matrix.txt"

if config["type"] == "single": # alignment
	rule fastq_to_sam:
		input:
			fastq = "fastq/{sample}_R1_001.fastq.gz"
		params:
			index = config["index"]
		output:
			sam = temp("processed/{sample}.sam")
		threads: config["alignment"]["threads"]
		log:
			"logs/{sample}.alignment.log"
		shell:
			"hisat2 -p {threads} -x {params.index} -U {input.fastq} "
			"-S {output.sam} 2> {log}"

elif config["type"] == "paired":
	rule fastq_to_sam:
		input:
			pair1 = "fastq/{sample}_R1_001.fastq.gz",
			pair2 = "fastq/{sample}_R2_001.fastq.gz"
		params:
			index = config["index"]
		output:
			sam = temp("processed/{sample}.sam")
		threads: config["alignment"]["threads"]
		log:
			"logs/{sample}.alignment.log"
		shell:
			"hisat2 -p {threads} -x {params.index} -1 {input.pair1} "
			"-2 {input.pair2} -S {output.sam} 2> {log}"

rule sam_to_unique: # chipseq
	input:
		sam = "processed/{sample}.sam"
	output:
		unique = temp("processed/{sample}.unique.sam")
	run:
		shell("grep -E 'NH:i:1|@' {input.sam} > {output.unique}")

if config["experiment"] == "rnaseq":
	rule sam_to_bam: # rnaseq
		input:
			sam = "processed/{sample}.sam"
		output:
			bam = temp("processed/{sample}.bam")
		shell:
			"samtools view -Sb {input.sam} > {output.bam}"

elif config["experiment"] == "chipseq":
	rule unique_to_bam: # chipseq
		input:
			sam = "processed/{sample}.unique.sam"
		output:
			bam = temp("processed/{sample}.bam")
		shell:
			"samtools view -Sb {input.sam} > {output.bam}"

rule bam_to_sortedbam:
	input:
		bam = "processed/{sample}.bam"
	output:
		sorted_bam = "processed/{sample}.sorted.bam"
	shell:
		"samtools sort -T processed/{wildcards.sample} "
		"-O bam {input.bam} > {output.sorted_bam}"

# chipseq 

rule sortedbam_to_rmdup:
	input:
		sorted_bam = "processed/{sample}.sorted.bam"
	output:
		dup_removed = "processed/{sample}.unique.sorted.rmdup.bam"
	log:
		"logs/{sample}.rmdup.log"
	run:
		shell("samtools rmdup {input.sorted_bam} {output.dup_removed} 2> {log}")
		shell("rm {input.sorted_bam}")

rule rmdup_to_tdf:
	input:
		dup_removed = "processed/{sample}.unique.sorted.rmdup.bam"
	params:
		chr_sizes = config["chr_sizes"]
	output:
		tdf = "processed/{sample}.unique.sorted.rmdup.tdf"
	log:
		"logs/{sample}.tdf.log"
	shell:
		"igvtools count {input.dup_removed} {output.tdf} {params.chr_sizes} "
		"2> {log}"

rule rmdup_to_chrbam:
	input:
		dup_removed = "processed/{sample}.unique.sorted.rmdup.bam"
	params:
		sam_chr_header = config["sam_chr_header"]
	output:
		chrbam = temp("processed/{sample}.unique.sorted.rmdup.chr.bam")
	log:
		"logs/{sample}.chrbam.log"
	shell:
		"samtools reheader {params.sam_chr_header} {input.dup_removed} > {output.chrbam}"

rule chrbam_to_bed:
	input:
		chrbam = "processed/{sample}.unique.sorted.rmdup.chr.bam"
	output:
		bed = "processed/{sample}.unique.sorted.rmdup.chr.bed"
	log:
		"logs/{sample}.bed.log"
	shell:
		"bedtools bamtobed -i {input.chrbam} > {output.bed} 2> {log}"

# rnaseq

rule sortedbam_to_counts:
	input:
		sorted_bam = "processed/{sample}.sorted.bam",
		gtf = config["gtf"]
	output:
		counts = "processed/{sample}.counts.txt"
	log:
		"logs/{sample}.htseq_counts.log"
	shell:
		"htseq-count -f bam -s no {input.sorted_bam} "
		"{input.gtf} > {output.counts} 2> {log}"

rule counts_matrix:
	input:
		counts = expand("processed/{sample}.counts.txt", \
						sample=SAMPLES)
	output:
		matrix = "processed/htseq_counts_matrix.txt"
	run:
		import pandas as pd

		dict_of_counts = {}

		for file in input:
			sample = file.split("/")[1].split(".")[0]
			dict_of_counts[sample] = {}

			with open(file, "r") as infile:
				for lines in infile:
					lines = lines.strip().split("\t")
					if "__" not in lines[0]:
						dict_of_counts[sample][lines[0]] = lines[1]

		dataframe = pd.DataFrame(dict_of_counts)
		dataframe.to_csv(output[0], sep = '\t')
