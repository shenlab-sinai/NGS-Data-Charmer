SAMPLES, = glob_wildcards("fastq/{sample}_R1_001.fastq.gz")
configfile: "config.yaml"

rule all:
	input:
		"counts/htseq_counts_matrix.txt"

if config["type"] == "single":
	rule alignment:
		input:
			fastq = "fastq/{sample}_R1_001.fastq.gz"
		params:
			index = config["index"]
		output:
			sam = temp("sorted_bam/{sample}.sam")
		threads: config["alignment"]["threads"]
		log:
			"logs/{sample}.alignment.log"
		shell:
			"hisat2 -p {threads} -x {params.index} -U {input.fastq} "
			"-S {output.sam} 2> {log}"

elif config["type"] == "paired":
	rule alignment:
		input:
			pair1 = "fastq/{sample}_R1_001.fastq.gz",
			pair2 = "fastq/{sample}_R2_001.fastq.gz"
		params:
			index = config["index"]
		output:
			sam = temp("sorted_bam/{sample}.sam")
		threads: config["alignment"]["threads"]
		log:
			"logs/{sample}.alignment.log"
		shell:
			"hisat2 -p {threads} -x {params.index} -1 {input.pair1} "
			"-2 {input.pair2} -S {output.sam} 2> {log}"

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
