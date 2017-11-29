configfile: "config.yaml"

rule all:
	input:
		"counts/htseq_counts_matrix.txt"

rule alignment:
	input:
		fastq = "fastq/{sample}.fastq"
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
		counts = temp("counts/{sample}.counts.txt")
	log:
		"logs/{sample}.htseq_counts.log"
	shell:
		"htseq-count -f bam -s no {input.sorted_bam} "
		"{input.gtf} > {output.counts} 2> {log}"

rule counts_matrix:
	input:
		counts = expand("counts/{sample}.counts.txt", \
						sample=config["samples"])
	output:
		matrix = "counts/htseq_counts_matrix.txt"
	run:
		import pandas as pd

		dict_of_counts = {}

		for file in input:
			sample = file.split("/")[1].split(".")[0]
			dict_of_counts[sample] = {}

			with open(file, "r") as infile, \
				open("logs/"+sample+".htseq_counts.log", "a") as outfile:
				for lines in infile:
					lines = lines.strip().split("\t")
					if "__" not in lines[0]:
						dict_of_counts[sample][lines[0]] = lines[1]
					else:
						outfile.write("\t".join(lines) + "\n")

		dataframe = pd.DataFrame(dict_of_counts)
		dataframe.to_csv(output[0], sep = '\t')
