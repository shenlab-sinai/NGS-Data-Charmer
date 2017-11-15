configfile: "config.yaml"

rule all:
	input:
		expand('counts/{sample}.counts.txt', sample=config["samples"])

rule alignment:
	input:
		fastq = "fastq/{sample}.fastq"
	params:
		index = config["index"]
	output:
		sam = temp("sorted_bam/{sample}.sam")
	threads: 2
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
		counts = "counts/{sample}.counts.txt"
	log:
		"logs/{sample}.htseq_counts.log"
	shell:
		"htseq-count -f bam -s no {input.sorted_bam} "
		"{input.gtf} > {output.counts} 2> {log}"
