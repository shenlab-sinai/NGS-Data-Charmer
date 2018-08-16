SAMPLES, = glob_wildcards("fastq/{sample}_R1_001.fastq.gz")
configfile: "config.yaml"

ALL_TDF = expand('processed/tdf/{sample}.unique.sorted.rmdup.tdf', sample=SAMPLES)
ALL_BW = expand('processed/bw/{sample}.unique.sorted.rmdup.chr.bw', sample=SAMPLES)
ALL_BED = expand('processed/bed/{sample}.unique.sorted.rmdup.chr.bed', sample=SAMPLES)

COUNTS_MATRIX = "processed/htseq_counts_matrix.txt"
MULTIQC_REPORT = "multiqc_report.html"

if config["experiment"] == "chipseq": # defining target
	rule all:
		input: ALL_BED, ALL_TDF, ALL_BW, MULTIQC_REPORT

elif config["experiment"] == "rnaseq":
	rule all:
		input: COUNTS_MATRIX, MULTIQC_REPORT

if config["type"] == "single": # alignment
	rule trim_fastq_fastqc:
		input:
			fastq = "fastq/{sample}_R1_001.fastq.gz"
		output:
			trimmed_fastq = temp("logs/{sample}_R1_001_trimmed.fq.gz"),
			fastqc_zipfile = "fastqc/{sample}_R1_001_fastqc.zip"
		log:
			"logs/{sample}.trim_adapters.log"
		run:
			shell("trim_galore {input.fastq} -o ./logs")
			shell("fastqc {input.fastq} -o ./fastqc")

	rule fastq_to_sam:
		input:
			trimmed_fastq = "logs/{sample}_R1_001_trimmed.fq.gz"
		params:
			index = config["index"]
		output:
			sam = temp("processed/bam/{sample}.sam")
		threads: config["threads_for_alignment"]
		log:
			"logs/{sample}.alignment.log"
		shell:
			"hisat2 -p {threads} -x {params.index} -U {input.trimmed_fastq} "
			"-S {output.sam} 2> {log}"

elif config["type"] == "paired":
	rule trim_fastq_fastqc:
		input:
			pair1 = "fastq/{sample}_R1_001.fastq.gz",
			pair2 = "fastq/{sample}_R2_001.fastq.gz"
		output:
			trimmed_pair1 = temp("logs/{sample}_R1_001_val_1.fq.gz"),
			trimmed_pair2 = temp("logs/{sample}_R2_001_val_2.fq.gz"),
			fastqc_zipfile1 = "fastqc/{sample}_R1_001_fastqc.zip",
			fastqc_zipfile2 = "fastqc/{sample}_R2_001_fastqc.zip"
		log:
			"logs/{sample}.trim_adapters.log"
		run:
			shell("trim_galore {input.pair1} {input.pair2} --paired -o ./logs")
			shell("fastqc {input.pair1} {input.pair2} -o ./fastqc")

	rule fastq_to_sam:
		input:
			trimmed_pair1 = "logs/{sample}_R1_001_val_1.fq.gz",
			trimmed_pair2 = "logs/{sample}_R2_001_val_2.fq.gz"
		params:
			index = config["index"]
		output:
			sam = temp("processed/bam/{sample}.sam")
		threads: config["threads_for_alignment"]
		log:
			"logs/{sample}.alignment.log"
		shell:
			"hisat2 -p {threads} -x {params.index} -1 {input.trimmed_pair1} "
			"-2 {input.trimmed_pair2} -S {output.sam} 2> {log}"

rule sam_to_unique: # chipseq
	input:
		sam = "processed/bam/{sample}.sam"
	output:
		unique = temp("processed/bam/{sample}.unique.sam")
	run:
		shell("grep -E 'NH:i:1|@' {input.sam} > {output.unique}")

rule sam_unique_to_bam:
	input:
		sam = "processed/bam/{sample}.sam" if config["experiment"] == "rnaseq" else "processed/bam/{sample}.unique.sam"
	output:
		bam = temp("processed/bam/{sample}.bam")
	shell:
		"samtools view -Sb {input.sam} > {output.bam}"

rule bam_to_sortedbam:
	input:
		bam = "processed/bam/{sample}.bam"
	output:
		sorted_bam = "processed/bam/{sample}.sorted.bam"
	shell:
		"samtools sort -T processed/bam/{wildcards.sample} "
		"-O bam {input.bam} > {output.sorted_bam}"

# chipseq 

rule sortedbam_to_rmdup:
	input:
		sorted_bam = "processed/bam/{sample}.sorted.bam"
	output:
		dup_removed = "processed/bam/{sample}.unique.sorted.rmdup.bam"
	log:
		"logs/{sample}.rmdup.log"
	run:
		shell("samtools rmdup {input.sorted_bam} {output.dup_removed} 2> {log}")
		shell("rm -f {input.sorted_bam}")

rule rmdup_to_tdf:
	input:
		dup_removed = "processed/bam/{sample}.unique.sorted.rmdup.bam"
	params:
		chr_sizes = config["chr_sizes"]
	output:
		tdf = "processed/tdf/{sample}.unique.sorted.rmdup.tdf"
	log:
		"logs/{sample}.tdf.log"
	shell:
		"igvtools count {input.dup_removed} {output.tdf} {params.chr_sizes} "
		"2> {log}"

rule rmdup_to_chrbam:
	input:
		dup_removed = "processed/bam/{sample}.unique.sorted.rmdup.bam"
	output:
		chrbam = "processed/bam/{sample}.unique.sorted.rmdup.chr.bam"
	log:
		"logs/{sample}.chrbam.log"
	run:
		shell('samtools view -H {input.dup_removed} | sed -e "s/SN:\([0-9XY]\)/SN:chr\\1/" -e "s/SN:MT/SN:chrM/" | samtools reheader - {input.dup_removed} > {output.chrbam}')

rule chrbam_to_bw:
	input:
		chrbam = "processed/bam/{sample}.unique.sorted.rmdup.chr.bam"
	output:
		bw_file = "processed/bw/{sample}.unique.sorted.rmdup.chr.bw"
	log:
		"logs/{sample}.bw.log"
	run:
		shell("samtools index {input.chrbam}")
		shell("bamCoverage -b {input.chrbam} -o {output.bw_file} --binSize 10 --normalizeUsing RPKM")

rule chrbam_to_bed:
	input:
		chrbam = "processed/bam/{sample}.unique.sorted.rmdup.chr.bam"
	output:
		bed = "processed/bed/{sample}.unique.sorted.rmdup.chr.bed"
	log:
		"logs/{sample}.bed.log"
	shell:
		"bedtools bamtobed -i {input.chrbam} > {output.bed} 2> {log}"

# rnaseq

rule sortedbam_to_counts:
	input:
		sorted_bam = "processed/bam/{sample}.sorted.bam",
		gtf = config["gtf"]
	output:
		counts = "processed/counts/{sample}.counts.txt"
	log:
		"logs/{sample}.htseq_counts.log"
	shell:
		"htseq-count -f bam -s no {input.sorted_bam} "
		"{input.gtf} > {output.counts} 2> {log}"

rule counts_matrix:
	input:
		counts = expand("processed/counts/{sample}.counts.txt", \
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

# multiqc
if config["experiment"] == "rnaseq":
	rule run_multiqc:
		input:
			matrix = "processed/htseq_counts_matrix.txt"
		output:
			multiqc_report = "multiqc_report.html"
		params:
			multiqc_config = config["multiqc_yaml"]
		shell:
			"multiqc . -f --config {params.multiqc_config}"

elif config["experiment"] == "chipseq":
	rule run_multiqc:
		input:
			tdf = expand("processed/tdf/{sample}.unique.sorted.rmdup.tdf", \
						sample = SAMPLES)
		output:
			multiqc_report = "multiqc_report.html"
		params:
			multiqc_config = config["multiqc_yaml"]
		shell:
			"multiqc . -f --config {params.multiqc_config}"
			
			
