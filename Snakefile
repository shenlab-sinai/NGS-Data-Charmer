SAMPLES, = glob_wildcards("fastq/{sample}_R1_001.fastq.gz")
configfile: "config.yaml"

def fastq_files(wildcards):
	if (config["type"] == "single"):
		return expand("fastq/{sample}_{strand}_001.fastq.gz", \
			strand=["R1"], sample=wildcards.sample)
	elif (config["type"] == "paired"):
		return expand("fastq/{sample}_{strand}_001.fastq.gz", \
			strand=["R1", "R2"], sample=wildcards.sample)

rule all:
	input:
		"counts/htseq_counts_matrix.txt"

rule alignment:
	input:
		fastq_files
	params:
		index = config["index"]
	output:
		sam = temp("sorted_bam/{sample}.sam")
	threads: config["alignment"]["threads"]
	log:
		"logs/{sample}.alignment.log"
	run:
		if config["type"] == "single":
			shell("hisat2 -p {threads} -x {params.index} -U \
				{input.sample} -S {output.sam} 2> {log}")
		elif config["type"] == "paired":
			shell("hisat2 -p {threads} -x {params.index} -1 \
				{input.forward} -2 {input.reverse} -S \
				{output.sam} 2> {log}")

rule sam_to_bam:
	input:
		sam = "sorted_bam/{sample}.sam"
	output:
		bam = temp("sorted_bam/{sample}.bam")
	shell:
		"samtools view -Sb {input.sam} > {output.bam}"

rule sort_bam:
	input:
		bam = "sorted_bam/{sample}.bam"
	output:
		sorted_bam = "sorted_bam/{sample}.sorted.bam"
	shell:
		"samtools sort -T sorted_bam/{wildcards.sample} "
		"-O bam {input.bam} > {output.sorted_bam}"

rule htseq_counts:
	input:
		sorted_bam = "sorted_bam/{sample}.sorted.bam",
		gtf = config["gtf"]
	output:
		counts = "counts/{sample}.counts.txt"
	log:
		"logs/{sample}.htseq_counts.log"
	shell:
		"htseq-count -f bam -s no {input.sorted_bam} "
		"{input.gtf} > {output.counts} 2> {log}"

rule counts_matrix:
	input:
		counts = expand("counts/{sample}.counts.txt", \
						sample=SAMPLES)
	output:
		matrix = "counts/htseq_counts_matrix.txt"
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
