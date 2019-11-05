SAMPLES, = glob_wildcards("fastq/{sample}_R1_001.fastq.gz")
configfile: "config.yaml"

ALL_TDF = expand('processed/tdf/{sample}.unique.sorted.rmdup.tdf', 
			sample=SAMPLES)
ALL_BW = expand('processed/bw/{sample}.unique.sorted.rmdup.chr.bw', 
			sample=SAMPLES)
ALL_BED = expand('processed/bed/{sample}.unique.sorted.rmdup.chr.bed', 
			sample=SAMPLES)

COUNTS_MATRIX = "processed/htseq_counts_matrix.txt"
MULTIQC_REPORT = "multiqc_report.html"

if config["experiment"] == "chipseq": # defining target
	rule all:
		input: ALL_BED, ALL_TDF, ALL_BW, MULTIQC_REPORT

elif config["experiment"] == "rnaseq":
	rule all:
		input: COUNTS_MATRIX, MULTIQC_REPORT
			
elif config["experiment"] == "cutrun":
    	rule all:
        	input: ALL_BED, ALL_BW, ALL_BAM, MULTIQC_REPORT

if config["experiment"] == "cutrun" and config["type"] == "paired":
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
    	rule split_length_long:
        	input:
            		trimg_pair1 = "logs/{sample}_R1_001_val_1.fq.gz",
            		trimg_pair2 = "logs/{sample}_R2_001_val_2.fq.gz"
        	output:
            		cut_r1_p1 = temp("logs/{sample}_t1_R1.len75.fastq"),
            		cut_r2_p1 = temp("logs/{sample}_t1_R2.len75.fastq")
        	params:
            		read_length = config["read_length"]
        	log:
            		"logs/{sample}.split_length_keeplong.log"
        	run:
            		shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r1_p1} {input.trimg_pair1}"),
            		shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r2_p1} {input.trimg_pair2}")

    	rule split_length_short:
        	input:
            		trimg_pair1 = "logs/{sample}_R1_001_val_1.fq.gz",
            		trimg_pair2 = "logs/{sample}_R2_001_val_2.fq.gz"
        	output:
            		cut_r1_p2 = temp("logs/{sample}_t1_R1.lt75.fastq"),
            		cut_r2_p2 = temp("logs/{sample}_t1_R2.lt75.fastq")
        	params:
            		read_length = config["read_length_max"]
        	log:
            		"logs/{sample}.split_length_keepshort.log"
        	run:
            		shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r1_p2} {input.trimg_pair1}"),
            		shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r2_p2} {input.trimg_pair2}")
    	rule trim_long:
        	input:
            		cut_r1_p1 = "logs/{sample}_t1_R1.len75.fastq",
            		cut_r2_p1 = "logs/{sample}_t1_R2.len75.fastq"
        	output:
            		cut_r1_p3 = temp("logs/{sample}_t1_R1.len75_trim.fastq"),
            		cut_r2_p3 = temp("logs/{sample}_t1_R2.len75_trim.fastq")
        	log:
            		"logs/{sample}.split_length_keepshort.log"
        	run:
            		shell("cutadapt -u -6 -o {output.cut_r1_p3} {input.cut_r1_p1}"),
            		shell("cutadapt -u -6 -o {output.cut_r2_p3} {input.cut_r2_p1}")

    	rule combine_split_lengths:
        	input:
            		cut_r1_p3 = "logs/{sample}_t1_R1.len75_trim.fastq",
            		cut_r2_p3 = "logs/{sample}_t1_R2.len75_trim.fastq",
            		cut_r1_p2 = "logs/{sample}_t1_R1.lt75.fastq",
            		cut_r2_p2 = "logs/{sample}_t1_R2.lt75.fastq"
        	output:
            		cut_r1_p4 = temp("logs/{sample}_t2_R1.fastq"),
            		cut_r2_p4 = temp("logs/{sample}_t2_R2.fastq")
        	run:
            		shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > {output.cut_r1_p4}"),
            		shell("cat {input.cut_r2_p3} {input.cut_r2_p2} > {output.cut_r2_p4}")

    	rule sort_combined_lengths:
        	input:
            		cut_r1_p4 = "logs/{sample}_t2_R1.fastq",
            		cut_r2_p4 = "logs/{sample}_t2_R2.fastq"
        	output:
            		cut_r1_p5 = temp("logs/{sample}_t2_R1_sorted.fastq.gz"),
            		cut_r2_p5 = temp("logs/{sample}_t2_R2_sorted.fastq.gz")
        	run:
            		shell("cat {input.cut_r1_p4} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\\n\" > logs/{wildcards.sample}_t2_R1_sorted.fastq"),
            		shell("cat {input.cut_r2_p4} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\\n\" > logs/{wildcards.sample}_t2_R2_sorted.fastq"),
            		shell("gzip logs/{wildcards.sample}_t2_R1_sorted.fastq"),
            		shell("gzip logs/{wildcards.sample}_t2_R2_sorted.fastq")

	rule bowtie2:
        	input:
            		trimmed_pair1 = "logs/{sample}_t2_R1_sorted.fastq.gz",
            		trimmed_pair2 = "logs/{sample}_t2_R2_sorted.fastq.gz"
        	params:
            		index = config["bowtie2_index"]
        	output:
            		bam = temp("processed/bam/{sample}.sorted.bam"),
            		bambai = temp("processed/bam/{sample}.sorted.bam.bai")
        	threads:
            		config["threads_for_alignment"]
        	log:
            		"logs/{sample}.alignment.log"
        	run:
            		shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} 2> {log} > processed/bam/{wildcards.sample}.sam"),
            		shell("samtools sort processed/bam/{wildcards.sample}.sam | samtools view -bS - > processed/bam/{wildcards.sample}.bam"),
            		shell("rm processed/bam/{wildcards.sample}.sam"),
            		shell("samtools view -bh -f 3 -F 4 -F 8 processed/bam/{wildcards.sample}.bam > processed/bam/{wildcards.sample}_mapped.bam"),
            		shell("samtools index processed/bam/{wildcards.sample}_mapped.bam"),
            		shell("samtools sort processed/bam/{wildcards.sample}_mapped.bam > {output.bam}"),
            		shell("samtools index {output.bam}")
            		shell("rm processed/bam/{wildcards.sample}_mapped.bam*")

	rule rmdup:
        	input:
            		"processed/bam/{sample}.sorted.bam"
        	output:
            		bam = "processed/bam/{sample}.unique.sorted.rmdup.bam"
        	params:
            		picardmetric = "logs/{sample}.markdups.metrics.txt"
        	run:
            		shell("picard MarkDuplicates INPUT={input} OUTPUT={output.bam} VALIDATION_STRINGENCY=SILENT METRICS_FILE={params.picardmetric}"),
            		shell("samtools index {output.bam}")

if config["experiment"] == "cutrun" and config["type"] == "single":
	rule trim_fastq_fastqc:
		input:
			pair1 = "fastq/{sample}_R1_001.fastq.gz"
		output:
		    	trimmed_pair1 = temp("logs/{sample}_R1_001_val_1.fq.gz"),
		    	fastqc_zipfile1 = "fastqc/{sample}_R1_001_fastqc.zip"
		log:
		    	"logs/{sample}.trim_adapters.log"
		run:
		    	shell("trim_galore {input.pair1} --paired -o ./logs")
		    	shell("fastqc {input.pair1} -o ./fastqc")
    	rule split_length_long:
        	input:
            		trimg_pair1 = "logs/{sample}_R1_001_val_1.fq.gz"
        	output:
            		cut_r1_p1 = temp("logs/{sample}_t1_R1.len75.fastq")
        	params:
            		read_length = config["read_length"]
        	log:
            		"logs/{sample}.split_length_keeplong.log"
        	run:
            		shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r1_p1} {input.trimg_pair1}")

    	rule split_length_short:
        	input:
            		trimg_pair1 = "logs/{sample}_R1_001_val_1.fq.gz"
        	output:
            		cut_r1_p2 = temp("logs/{sample}_t1_R1.lt75.fastq")
        	params:
            		read_length = config["read_length_max"]
        	log:
            		"logs/{sample}.split_length_keepshort.log"
        	run:
            		shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r1_p2} {input.trimg_pair1}")
    	rule trim_long:
        	input:
            		cut_r1_p1 = "logs/{sample}_t1_R1.len75.fastq"
        	output:
            		cut_r1_p3 = temp("logs/{sample}_t1_R1.len75_trim.fastq")
        	log:
            		"logs/{sample}.split_length_keepshort.log"
        	run:
            		shell("cutadapt -u -6 -o {output.cut_r1_p3} {input.cut_r1_p1}")

    	rule combine_split_lengths:
        	input:
            		cut_r1_p3 = "logs/{sample}_t1_R1.len75_trim.fastq",
            		cut_r1_p2 = "logs/{sample}_t1_R1.lt75.fastq"
        	output:
            		cut_r1_p4 = temp("logs/{sample}_t2_R1.fastq")
        	run:
            		shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > {output.cut_r1_p4}")			

	rule bowtie2:
        	input:
            		trimmed_pair1 = "logs/{sample}_t2_R1.fastq.gz"
        	params:
            		index = config["bowtie2_index"]
        	output:
            		bam = temp("processed/bam/{sample}.sorted.bam"),
            		bambai = temp("processed/bam/{sample}.sorted.bam.bai")
        	threads:
            		config["threads_for_alignment"]
        	log:
            		"logs/{sample}.alignment.log"
        	run:
            		shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -U {input.trimmed_pair1} 2> {log} > processed/bam/{wildcards.sample}.sam"),
            		shell("samtools sort processed/bam/{wildcards.sample}.sam | samtools view -bS - > processed/bam/{wildcards.sample}.bam"),
            		shell("rm processed/bam/{wildcards.sample}.sam"),
            		shell("samtools view -bh -f 3 -F 4 -F 8 processed/bam/{wildcards.sample}.bam > processed/bam/{wildcards.sample}_mapped.bam"),
            		shell("samtools index processed/bam/{wildcards.sample}_mapped.bam"),
            		shell("samtools sort processed/bam/{wildcards.sample}_mapped.bam > {output.bam}"),
            		shell("samtools index {output.bam}")
            		shell("rm processed/bam/{wildcards.sample}_mapped.bam*")

	rule rmdup:
        	input:
            		"processed/bam/{sample}.sorted.bam"
        	output:
            		bam = "processed/bam/{sample}.unique.sorted.rmdup.bam"
        	params:
            		picardmetric = "logs/{sample}.markdups.metrics.txt"
        	run:
            		shell("picard MarkDuplicates INPUT={input} OUTPUT={output.bam} VALIDATION_STRINGENCY=SILENT METRICS_FILE={params.picardmetric}"),
            		shell("samtools index {output.bam}")			
			
	
	
	
	
	
			
if config["type"] == "single" and config["experiment"] != "cutrun" : # alignment
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
			"hisat2 -p {threads} -x {params.index} -U "
			"{input.trimmed_fastq} -S {output.sam} 2> {log}"

elif config["type"] == "paired" and config["experiment"] != "cutrun" :
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
			"hisat2 -p {threads} -x {params.index} -1 "
			"{input.trimmed_pair1} -2 {input.trimmed_pair2} "
			"-S {output.sam} 2> {log}"
			
if config["experiment"] != "cutrun" :
	rule sam_to_unique: # chipseq
		input:
			sam = "processed/bam/{sample}.sam"
		output:
			unique = temp("processed/bam/{sample}.unique.sam")
		run:
			shell("grep -E 'NH:i:1|@' {input.sam} > {output.unique}")

	rule sam_unique_to_bam:
		input:
			sam = "processed/bam/{sample}.sam" if config["experiment"] == "rnaseq" 
					else "processed/bam/{sample}.unique.sam"
		output:
			bam = temp("processed/bam/{sample}.bam")
		shell:
			"samtools view -Sb {input.sam} > {output.bam}"

	rule bam_to_sortedbam:
		input:
			bam = "processed/bam/{sample}.bam"
		output:
			sorted_bam = "processed/bam/{sample}.sorted.bam"
		threads: config["threads_for_alignment"]
		shell:
			"samtools sort -T processed/bam/{wildcards.sample} "
			"-O bam -@ {threads} -o {output.sorted_bam} {input.bam}"

	if config["experiment"] == "rnaseq":
		rule sortedbam_to_rmdup_rnaseq:
			input:
				sorted_bam = "processed/bam/{sample}.sorted.bam"
			output:
				dup_removed = "processed/bam/{sample}.sorted.rmdup.bam"
			log:
				"logs/{sample}.rmdup.log"
			run:
				shell("samtools rmdup {input.sorted_bam} {output.dup_removed} 2> {log}")
				shell("rm -f {input.sorted_bam}")

	# chipseq 
	if config["experiment"] == "chipseq":
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
	shell:
		'samtools view -H {input.dup_removed} '
		'| sed -e "s/SN:\([0-9XY]\)/SN:chr\\1/" -e "s/SN:MT/SN:chrM/" '
		'| samtools reheader - {input.dup_removed} > {output.chrbam}'

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
if config["count_scheme"] == "fraction":
	rule sortedbam_to_counts:
		input:
			sorted_bam = "processed/bam/{sample}.sorted.bam" if config["type"] == "single"
							else "processed/bam/{sample}.sorted.rmdup.bam",
			gtf = config["gtf"]
		output:
			counts = "processed/counts/{sample}.counts.txt"
		log:
			"logs/{sample}.htseq_counts.log"
		shell:
			"featureCounts -p -O --fraction -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}" if config["type"] == "paired" 
			else "featureCounts -O --fraction -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}"

elif config["count_scheme"] == "count_all":
	rule sortedbam_to_counts:
		input:
			sorted_bam = "processed/bam/{sample}.sorted.bam" if config["type"] == "single"
							else "processed/bam/{sample}.sorted.rmdup.bam",
			gtf = config["gtf"]
		output:
			counts = "processed/counts/{sample}.counts.txt"
		log:
			"logs/{sample}.htseq_counts.log"
		shell:
			"featureCounts -p -O -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}" if config["type"] == "paired" 
			else "featureCounts -O -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}"

elif config["count_scheme"] == "count_uniq":
	rule sortedbam_to_counts:
		input:
			sorted_bam = "processed/bam/{sample}.sorted.bam" if config["type"] == "single"
							else "processed/bam/{sample}.sorted.rmdup.bam",
			gtf = config["gtf"]
		output:
			counts = "processed/counts/{sample}.counts.txt"
		log:
			"logs/{sample}.htseq_counts.log"
		shell:
			"featureCounts -p -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}" if config["type"] == "paired" 
			else "featureCounts -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}"

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
			sample = file.split(".")[0]
			dict_of_counts[sample] = {}

			with open(file, "r") as infile:
				next(infile)
				next(infile)
				for lines in infile:
					lines = lines.strip().split("\t")
					dict_of_counts[sample][lines[0]] = int(float(lines[6]))

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
			
			
