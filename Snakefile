import os
from os import listdir
from os.path import isfile, join
import re
from Bio import SeqIO
import gzip

configfile: "config.yaml"
myfastqpath = "fastq/"

# HELPER FUNCTIONS

# Create the pattern of file endings

def create_endings(x):
	"""
	Returns a list of likely fastq file endings

	Input Parameter: 
	x (int): Either 1 or 2, indicating the forward (1) or reverse (2) read.

	Returns: 
	list: A list of strings, representing the file endings a user might use for denoting their fastq files.
	"""
	return(["_R" + str(x) + "_001.fastq", "_R" + str(x) + "_001.fq", "_R" + str(x) + ".fastq", "_R" + str(x) + ".fq", "_" + str(x) + ".fastq", "_" + str(x) + ".fq", ".R" + str(x) + "_001.fastq", ".R" + str(x) + "_001.fq", ".R" + str(x) + ".fastq", ".R" + str(x) + ".fq", "." + str(x) + ".fastq", "." + str(x) + ".fq", "_r" + str(x) + "_001.fastq", "_r" + str(x) + "_001.fq", "_r" + str(x) + ".fastq", "_r" + str(x) + ".fq", ".r" + str(x) + "_001.fastq", ".r" + str(x) + "_001.fq", ".r" + str(x) + ".fastq", ".r" + str(x) + ".fq"])

# Function to list the fastq files present in the fastq folder
def getfilelist(myfastqpath):
	"""
	Extracts fastq files from the list of files present in your fastq directory.

	Input Parameter: 
	myfastqpath (string): directory containing your fastq files.
	
	Returns:
	list: List containing two strings. 
	1st string is all non-metadata files in the fastq directory
	2nd string is all non-metadata files ending in '.gz'
	"""
	onlyfiles = [f for f in listdir(myfastqpath) if isfile(join(myfastqpath, f))]
	onlyfiles = [i for i in onlyfiles if i.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))]
	gzfiles = [i for i in onlyfiles if i.endswith((".gz"))]
	return([onlyfiles, gzfiles])

# Unify fastq files to single file ending
def fix_input_files(file_suffix, input_fileset):
	"""
	Renames mixed input fastq files to the most common file ending and returns the selected file ending. NOTE: This step permenantly renames your fastq files from their original file ending. 

	Input Parameter: 
	file_suffix (string): ".gz" or ""; Gzipped fastq files are expected to end with the suffix ".gz". If files are NOT gzipped, the input is "".
	
	input_fileset (list): List of fastq file names to be examined. As written, gzipped files are listed within the variable 'gzfiles' and non-gzipped files are listed within the variable 'onlyfiles'.

	Returns: 
	list: A list containing four strings, the selected Read1 (forward read) file ending and the corresponding Read2 (reverse read) file ending, a list of all fastq-like files, and a list of gzipped fastq-like files.
	"""
	# Create the series of fastq file endings to search
	base_endings_r1 = create_endings(1)
	base_endings_r2 = create_endings(2)
	ending_dictionary = dict(zip(base_endings_r1, base_endings_r2))  # Define the R1 and R1 suffix pairs for reference
	mylist = list()  # Create empty list 
	for x in base_endings_r1:  # first traverse the R1 base endings to find the common ending
		matched_ends = [i for i in input_fileset if i.endswith(x + file_suffix)]
		if(len(matched_ends) > 0):
			mylist.append((x, len(matched_ends)))

	if len(mylist) == 0:  # Contingency: if all samples are single-end and do not have a common forward read ending
		print("Your dataset appears to be entirely single-end files.")
		odd_files = [i for i in input_fileset if i.endswith(".fq" + file_suffix)]
		if len(odd_files) > 0:
			oldnames = odd_files
			old_rep = [i.replace(".fq" + suffix, ".fastq" + suffix) for i in oldnames]
			[os.rename(join(myfastqpath, i),join(myfastqpath, y)) for i, y in zip(oldnames, old_rep)]
		# Re-assess and return filenames
		return([".fastq", ".fastq", getfilelist(myfastqpath)[0], getfilelist(myfastqpath)[1]])

	else:  # If R1 endings are present, check values and correct file names
		mylist_dict = {key: value for (key, value) in mylist}  # create dictionary of mixed file endings
		myR1_suffix = [key for (key, value) in mylist_dict.items() if value == max(mylist_dict.values())][0]
		myR2_suffix = ending_dictionary[myR1_suffix]
		mylist_dict.pop(myR1_suffix)  # remove main R1 suffix from dictionary
		if len(mylist_dict) > 0:  # begin processing forward reads
			for x in mylist_dict:
				oldnames = [i for i in input_fileset if i.endswith(x + file_suffix)]
				old_rep = [i.replace(x, myR1_suffix) for i in oldnames]
				[os.rename(join(myfastqpath, i),join(myfastqpath, y)) for i, y in zip(oldnames, old_rep)]

		mylist = list()  # Create empty list 
		for x in base_endings_r2:  # first traverse the R2 base endings to find the common ending
			matched_ends = [i for i in input_fileset if i.endswith(x + file_suffix)]
			if(len(matched_ends) > 0):
				mylist.append(x)
		mylist.remove(myR2_suffix)  # remove main R2 suffix from list
		if len(mylist) > 0:  # begin processing forward reads
			for x in mylist:
				oldnames = [i for i in input_fileset if i.endswith(x + file_suffix)]
				old_rep = [i.replace(x, myR2_suffix) for i in oldnames]
				[os.rename(join(myfastqpath, i),join(myfastqpath, y)) for i, y in zip(oldnames, old_rep)]

		#Re-assess file names
		if file_suffix == ".gz":
			input_fileset = getfilelist(myfastqpath)[1]
		else:
			input_fileset = getfilelist(myfastqpath)[0]

		# Now process single end files
		odd_files = [i for i in [i for i in input_fileset if not i.endswith(myR1_suffix + file_suffix)] if not i.endswith(myR2_suffix + file_suffix)]
		fastq_odd_1 = [i for i in odd_files if i.endswith(".fastq" + file_suffix)]
		fastq_odd_2 = [i for i in odd_files if i.endswith(".fq" + file_suffix)]
		if len(odd_files) > 0:
			print("Now unifying " + str(len(odd_files)) + " single-end files to \""+ myR1_suffix + file_suffix + "\" ending")
			if len(fastq_odd_1) > 0:  # rename files to correct ending
				oldnames = odd_files
				old_rep = [i.replace(".fastq" + file_suffix, myR1_suffix + file_suffix) for i in oldnames]
				[os.rename(join(myfastqpath, i),join(myfastqpath, y)) for i, y in zip(oldnames, old_rep)]
			if len(fastq_odd_2) > 0:  # rename files to correct ending
				oldnames = odd_files
				old_rep = [i.replace(".fq" + file_suffix, myR1_suffix + file_suffix) for i in oldnames]
				[os.rename(join(myfastqpath, i),join(myfastqpath, y)) for i, y in zip(oldnames, old_rep)]
		# Re-assess and return filenames and file endings
		return([myR1_suffix, myR2_suffix, getfilelist(myfastqpath)[0], getfilelist(myfastqpath)[1]])

# Function to retrieve and check cut&run read lengths
def check_readlength(suffix, input_fileset):
	"""
	When samples are specified to be cut&run:
	Samples the first read in each forward read file and extracts the read length. 
	Checks if the samples have different read lengths.

	Input Parameter: 
	suffix (string): ".gz" or ""; Gzipped fastq files are expected to end with the suffix ".gz". If files are NOT gzipped, the input is "".
	
	input_fileset (list): List of fastq file names to be examined. As written, gzipped files are listed within the variable 'gzfiles' and non-gzipped files are listed within the variable 'onlyfiles'.

	
	Returns:
	list: List containing two integers. 
	1st integer is the read length
	2nd integer is the read length, minus one
	"""
	my_cr_files = [i for i in input_fileset if i.endswith(R1_file_ending + suffix)]
	if suffix == ".gz":
		try:  # Extract the read length of first forward read
			read_len_list = [len(next(SeqIO.parse(gzip.open(join(myfastqpath, i), "rt"), "fastq")).seq) for i in my_cr_files]
			dedup_lengths = list(dict.fromkeys(read_len_list))
		except:
			raise NameError("One of your fastq files may be empty\nNow aborting...")
	elif suffix == "": 
		try:
			read_len_list = [len(next(SeqIO.parse(join(myfastqpath, i), "fastq")).seq) for i in my_cr_files]
			dedup_lengths = list(dict.fromkeys(read_len_list))
		except:
			raise NameError("One of your fastq files may be empty\nNow aborting...")

	if len(dedup_lengths) > 1:
		raise NameError("Based on sampling the first read of each R1 fastq file, your cut&run files have different read lengths!\nRecorded lengths:"+ ("elements in the list are "+', '.join(['%.f']*len(dedup_lengths))) % tuple(dedup_lengths) + " base pairs" + "\nAborting...")
	else:
		print("Congratulations, your cut&run fastq files appear to have uniform sequence lengths!\nProceeding with a read length of " + format(dedup_lengths[0]))
		return(dedup_lengths[0], dedup_lengths[0]-1)

# Create function for creating rule sets
def choose_rule_all(config):
	"""
	Selects the input needed for 'rule all' in snakemake pipeline

	Input Parameter: 
	config (dict): Dictionary derived from config.yaml and any additional key:value pairs added during the file preperation steps. 

	Returns:
	list: List of required input files for 'rule all'. 
	"""
	myout = []
	if config["to_multiqc"] == "TRUE":
		myout.append("output/multiqc_report.html")
	if config["to_bw"] == "TRUE" and config["experiment"] != "rnaseq":
		myout.append(expand('output/bw/{sample}.unique.sorted.rmdup.chr.bw', sample=SAMPLES))
	if config["to_bed"] == "TRUE" and config["experiment"] != "rnaseq":
		myout.append(expand('output/bed/{sample}.unique.sorted.rmdup.chr.bed', sample=SAMPLES))
	if config["to_tdf"] == "TRUE" and config["experiment"] != "rnaseq":
		myout.append(expand('output/tdf/{sample}.unique.sorted.rmdup.tdf', sample=SAMPLES))
	if config["experiment"] == "rnaseq":
		myout.append("output/htseq_counts_matrix.txt")
	# myout.append(expand('output/trim_fastq/{sample}_R1_001_trimmed.fq.gz', sample=SAMPLES))
	return(myout)

# Create read-pair inputs for sample processing
def create_inputs(config) :
	"""
	Creates the fastq file inputs needed for read trimming steps of the snakemake pipeline

	Input Parameter: 
	config (dict): Dictionary derived from config.yaml and any additional key:value pairs added during the file preperation steps. 

	Returns:
	list: List of two strings; 
	1st string denotes the forward read
	2nd string denotes the reverse read
	"""
	return([("fastq/{sample}" + expand("{ending}{suffix}",ending=R1_file_ending,suffix=suffix)[0]+""),("fastq/{sample}" + expand("{ending}{suffix}",ending=R2_file_ending,suffix=suffix)[0]+"")])

# END HELPER FUNCTIONS


# Retrieve the list of fastq files
onlyfiles = getfilelist(myfastqpath)[0]
gzfiles = getfilelist(myfastqpath)[1]

# Raise exception if no fastq files present
if len(gzfiles) == 0 and len(onlyfiles) == 0:
	raise NameError("You do not seem to have any fastq files present to process. Now exiting...")

# Raise exception if fastq files are a mixture of gzipped and non-gzipped files
if len(gzfiles) > 0 and len(gzfiles) != len(onlyfiles):
	myinput = "You have a mixture of gzipped files and non-gzipped files\nOnly {} of total {} files are gzipped!"
	raise NameError(print(myinput.format(len(gzfiles), len(onlyfiles))))

# Unify fastq file endings and return the final ending to be used.
if len(gzfiles) > 0:
	R1_file_ending, R2_file_ending, onlyfiles, gzfiles, = fix_input_files(".gz", gzfiles)
	suffix = ".gz"
else:
	R1_file_ending, R2_file_ending, onlyfiles, gzfiles = fix_input_files("", onlyfiles)
	suffix = ""

sample_string = myfastqpath + "{sample}" + R1_file_ending + suffix
SAMPLES, = glob_wildcards(sample_string)

# Check the file pairing
# Raise exception for non-paired PE files
if config["type"] == "single":
	print("You have chosen to use single-end reads\nRead pairing not being checked...")
elif config["type"] == "paired":
	len_r1 = len([i for i in onlyfiles if i.endswith(R1_file_ending + suffix)])
	if len_r1*2 != len(onlyfiles):
		myinput = "One or more samples do not have a read pair!\nIf using paired-end samples, please ensure each sample has read 1 and read 2 files\nAborting..."
		raise NameError(myinput) # Raise exception to break workflow
else:
	myinput = "You have specified unknown read type: " + config["type"] + "\nPlease specify either \"paired\" or \"single\" in the config.yaml file, then rerun the pipeline."
	raise NameError(myinput)

# Generate input rule for Snakemake
print(choose_rule_all(config))

rule all :
	input:  
		choose_rule_all(config)

# Begin Snakemake pre-processing for Cut&Run samples

subworkflow cutrun_preprocess:
	snakefile:
		"Snakefile_CR_preprocessing"
	configfile:
		"config.yaml"

if config["experiment"] == "cutrun" :
	rule bowtie2:
		input:
			trimmed_pair1 = cutrun_preprocess("output/trim_fastq/{sample}_R1_trimmed.fq.gz"),
			trimmed_pair2 = cutrun_preprocess("output/trim_fastq/{sample}_R2_trimmed.fq.gz")
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
			if config["type"] == "paired" :
				shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} 2> {log} > output/bam/{wildcards.sample}.sam"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > output/bam/{wildcards.sample}.bam"),
				shell("rm output/bam/{wildcards.sample}.sam"),
				shell("samtools index output/bam/{wildcards.sample}.bam")
			else : 
				shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -U {input.trimmed_pair1} 2> {log} > output/bam/{wildcards.sample}.sam"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > output/bam/{wildcards.sample}.bam"),
				shell("rm output/bam/{wildcards.sample}.sam"),
				shell("samtools index output/bam/{wildcards.sample}.bam")
				shell("rm {input.trimmed_pair2}")
			if config["keep_fastq"] == "FALSE" and config["type"] == "paired":
				shell("rm {input.trimmed_pair1} {input.trimmed_pair2}")
			elif config["keep_fastq"] == "FALSE" and config["type"] == "single":
				shell("rm {input.trimmed_pair1}")
			# Now sorting the bam file
			shell("samtools view -bh -f 3 -F 4 -F 8 output/bam/{wildcards.sample}.bam > output/bam/{wildcards.sample}_mapped.bam")
			shell("samtools index output/bam/{wildcards.sample}_mapped.bam")
			shell("samtools sort output/bam/{wildcards.sample}_mapped.bam > output/bam/{wildcards.sample}.sorted.bam")
			shell("samtools index output/bam/{wildcards.sample}.sorted.bam")
			shell("rm output/bam/{wildcards.sample}_mapped.bam*")
			if config["keep_unfiltered_bam"] == "FALSE":
				shell("rm output/bam/{wildcards.sample}.bam")
				shell("rm output/bam/{wildcards.sample}.bam.bai")

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


if config["experiment"] == "rnaseq" or config["experiment"] == "chipseq":
	subworkflow RNAandCHIP_preprocess:
		snakefile:
			"Snakefile_RNAandCHIP_preprocessing"
		configfile:
			"config.yaml"

	if config["experiment"] == "chipseq":
		rule bam_to_unique_mapped:	##Remove the multimapped reads and sort
			input:
				RNAandCHIP_preprocess("output/bam/{sample}.bam")
			output:
				bam = temp("output/bam/{sample}.sorted.bam")
			run:
				if config["type"] == "paired" :
					shell("samtools view -bh -f 3 -F 4 -F 8 -F 256 {input} > output/bam/{wildcards.sample}_filtered.bam") # mapped with pair, unique
					shell("samtools sort -O BAM -o {output.bam} output/bam/{wildcards.sample}_filtered.bam")
					shell("rm output/bam/{wildcards.sample}_filtered.bam")
				else : 
					shell("samtools view -bh -F 4 -F 256 {input} > output/bam/{wildcards.sample}_filtered.bam") # mapped with pair, unique locations
					shell("samtools sort -O BAM -o {output.bam} output/bam/{wildcards.sample}_filtered.bam")
					shell("rm output/bam/{wildcards.sample}_filtered.bam")
				if config["keep_unfiltered_bam"] == "FALSE" :
					shell("rm output/bam/{wildcards.sample}.bam output/bam/{wildcards.sample}.bam.bai")

		rule sortedbam_to_rmdup:
			input:
				sorted_bam = "output/bam/{sample}.sorted.bam"
			output:
				dup_removed = "output/bam/{sample}.unique.sorted.rmdup.bam"
			log:
				"output/logs/{sample}.rmdup.log"
			run:
				shell("samtools rmdup {input.sorted_bam} {output.dup_removed} 2> {log}")
				if config["keep_unfiltered_bam"] == "FALSE" :
					shell("rm -f {input.sorted_bam} {input.sorted_bam}.bai")			

	elif config["experiment"] == "rnaseq" and config["type"] == "paired" :
		rule sortedbam_to_rmdup:
			input:
				RNAandCHIP_preprocess("output/bam/{sample}.bam")
			output:
				"output/bam/{sample}.sorted.rmdup.bam"
			log:
				"output/logs/{sample}.rmdup.log"
			run:
				shell("samtools rmdup {input} {output} 2> {log}")
				if config["keep_unfiltered_bam"] == "FALSE" :
					shell("rm -f {input} {input}.bai")

	elif config["experiment"] == "rnaseq" and config["type"] == "single" :
		rule sortedbam_to_rmdup:
			input:
				RNAandCHIP_preprocess("output/bam/{sample}.bam")
			output:
				"output/bam/{sample}.sorted.bam"
			log:
				"output/logs/{sample}.rmdup.log"
			run:
				shell("cp {input} {output}")
				if config["keep_unfiltered_bam"] == "FALSE" :
					shell("rm -f {input} {input}.bai")

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

rule rmdup_to_chrbam:
	input:
		dup_removed = "output/bam/{sample}.unique.sorted.rmdup.bam"
	output:
		chrbam = "output/bam/{sample}.unique.sorted.rmdup.chr.bam",
		chrbambai = "output/bam/{sample}.unique.sorted.rmdup.chr.bam.bai"
	log:
		"output/logs/{sample}.chrbam.log"
	run:
		shell("samtools view -H {input.dup_removed} | sed -e \"s/SN:\([0-9XY]\)/SN:chr\\1/\" -e \"s/SN:MT/SN:chrM/\" | samtools reheader - {input.dup_removed} > {output.chrbam}")
		shell("samtools index {output.chrbam}")

rule chrbam_to_bw:
	input:
		chrbam = "output/bam/{sample}.unique.sorted.rmdup.chr.bam",
		chrbambai = "output/bam/{sample}.unique.sorted.rmdup.chr.bam.bai"
	output:
		bw_file = "output/bw/{sample}.unique.sorted.rmdup.chr.bw"
	log:
		"output/logs/{sample}.bw.log"
	run:
		shell("bamCoverage -b {input.chrbam} -o {output.bw_file} --binSize 10 --normalizeUsing RPKM")

rule chrbam_to_bed:
	input:
		chrbam = "output/bam/{sample}.unique.sorted.rmdup.chr.bam"
	output:
		bed = "output/bed/{sample}.unique.sorted.rmdup.chr.bed"
	log:
		"output/logs/{sample}.bed.log"
	shell:
		"bedtools bamtobed -i {input.chrbam} > {output.bed} 2> {log}"
							
# rnaseq
config["gene_scheme"] = "-t gene" ## or "-g gene_id"
rule sortedbam_to_counts:
	input:
		sorted_bam = "output/bam/{sample}.sorted.bam" if config["type"] == "single" else "output/bam/{sample}.sorted.rmdup.bam"
	output:
		counts = "output/counts/{sample}.counts.txt"
	params :
		gtf = config["gtf"],
		gene_scheme = config["gene_scheme"]
	log:
		"output/logs/{sample}.htseq_counts.log"
	run :
		if config["count_scheme"] == "fraction" and config["type"] == "paired" :
			shell("featureCounts -p -O --fraction {params.gene_scheme} -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}") 
		elif config["count_scheme"] == "fraction" and config["type"] == "single" :
			shell("featureCounts -O --fraction {params.gene_scheme} -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}") 
		elif config["count_scheme"] == "count_all" and config["type"] == "paired" :
			shell("featureCounts -p -O {params.gene_scheme} -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}") 
		elif config["count_scheme"] == "count_all" and config["type"] == "single" :
			shell("featureCounts -O {params.gene_scheme} -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}") 
		elif config["count_scheme"] == "count_uniq" and config["type"] == "paired" :
			shell("featureCounts -p {params.gene_scheme} -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}") 
		elif config["count_scheme"] == "count_uniq" and config["type"] == "single" :
			shell("featureCounts {params.gene_scheme} -a {params.gtf} -o {output.counts} {input.sorted_bam} 2> {log}")

rule counts_matrix:
	input:
		counts = expand("output/counts/{sample}.counts.txt", sample=SAMPLES)
	output:
		matrix = "output/htseq_counts_matrix.txt"
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
					if {params} == "-t gene" :
						dict_of_counts[sample][lines[0]] = int(float(lines[7]))
					else :
						dict_of_counts[sample][lines[0]] = int(float(lines[6]))

		dataframe = pd.DataFrame(dict_of_counts)
		dataframe.to_csv(output[0], sep = '\t')			

# multiqc
if config["experiment"] == "rnaseq":
	rule run_multiqc:
		input:
			matrix = "output/htseq_counts_matrix.txt"
		output:
			multiqc_report = "output/multiqc_report.html"
		params:
			multiqc_config = config["multiqc_yaml"]
		shell:
			"multiqc . -f --config {params.multiqc_config}"

elif config["experiment"] == "chipseq" or config["experiment"] == "cutrun" :
	rule run_multiqc:
		input:
			chrbam = expand("output/bam/{sample}.unique.sorted.rmdup.chr.bam", sample = SAMPLES)
		output:
			multiqc_report = "output/multiqc_report.html"
		params:
			multiqc_config = config["multiqc_yaml"]
		shell:
			"multiqc . -f --config {params.multiqc_config}"
