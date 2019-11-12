configfile: "config.yaml"

myfastqpath="./fastq/"
import os
from os import listdir
from os.path import isfile, join
import re
from Bio import SeqIO
import gzip
onlyfiles = [f for f in listdir(myfastqpath) if isfile(join(myfastqpath, f))]
gzfiles = list(filter(lambda x: re.search(r'.gz$', x), onlyfiles))

if len(gzfiles)>0 :
	if len(gzfiles)!=len(onlyfiles) : 	##Error message:
		myinput="You have a mixture of gzipped files and non-gzipped files\nOnly {} of total {} files are gzipped!"
		raise NameError(print(myinput.format(len(gzfiles),len(onlyfiles)))) ##exception to break workflow

base_endings_r1 = ["_R1_001.fastq", "_R1_001.fq", "_R1.fastq", "_R1.fq", "_1.fastq","_1.fq", ".R1_001.fastq", ".R1_001.fq", ".R1.fastq", ".R1.fq", ".1.fastq", ".1.fq","_r1_001.fastq", "_r1_001.fq", "_r1.fastq", "_r1.fq", ".r1_001.fastq", ".r1_001.fq", ".r1.fastq", ".r1.fq"]
base_endings_r2 = ["_R2_001.fastq", "_R2_001.fq", "_R2.fastq", "_R2.fq", "_2.fastq","_2.fq", ".R2_001.fastq", ".R2_001.fq", ".R2.fastq", ".R2.fq", ".2.fastq", ".2.fq","_r2_001.fastq", "_r2_001.fq", "_r2.fastq", "_r2.fq", ".r2_001.fastq", ".r2_001.fq", ".r2.fastq", ".r2.fq"]


def fix_input_files_R1(base_endings,file_suffix,input_fileset,return_name) :
	mylist = list()
	for x in base_endings :
		my_regex = re.escape(x) + file_suffix
		if(len([i for i in input_fileset if re.search(my_regex, i)])>0) :
			mylist.append((x,len([i for i in input_fileset if re.search(my_regex, i)])))
	mylist_dict = {key: value for (key, value) in mylist}
	inverse = [(value, key) for key, value in mylist_dict.items()]
	suffixes_to_change = [x for x in list(mylist_dict.keys()) if x!=max(inverse)[1]]
	if len(suffixes_to_change)>0 :
		for x in suffixes_to_change :
			my_regex = re.escape(x) + file_suffix
			oldnames = [i for i in input_fileset if re.search(my_regex, i)]
			print(oldnames)
			for myseq in range(len(oldnames)) :
				print(oldnames[myseq])
				print([i.replace(x,max(inverse)[1]) for i in oldnames][myseq])
				os.rename(join(myfastqpath,oldnames[myseq]),join(myfastqpath,[i.replace(x,max(inverse)[1]) for i in oldnames][myseq]))
	if return_name == "TRUE" :
		return max(inverse)[1]


def fix_input_files_R2(base_endings,file_suffix,input_fileset,r2_name) :
	mylist = list()
	for x in base_endings :
		my_regex = re.escape(x) + file_suffix
		if(len([i for i in input_fileset if re.search(my_regex, i)])>0) :
			mylist.append((x,len([i for i in input_fileset if re.search(my_regex, i)])))
	mylist_dict = dict(mylist)
	suffixes_to_change = [x for x in list(mylist_dict.keys()) if x!=r2_name]
	if len(suffixes_to_change)>0 :
		for x in suffixes_to_change :
			my_regex = re.escape(x) + file_suffix
			oldnames = [i for i in input_fileset if re.search(my_regex, i)]
			print(oldnames)
			for myseq in range(len(oldnames)) :
				print(oldnames[myseq])
				print([i.replace(x,r2_name) for i in oldnames][myseq])
				os.rename(join(myfastqpath,oldnames[myseq]),join(myfastqpath,[i.replace(x,r2_name) for i in oldnames][myseq]))

if len(gzfiles)>0 :
	base_endings_r1_complete = base_endings_r1 + [s + ".gz" for s in base_endings_r1]
	base_endings_r2_complete = base_endings_r2 + [s + ".gz" for s in base_endings_r2]
	R1_file_ending = fix_input_files_R1(base_endings_r1,".gz$",gzfiles,"TRUE")
	R2_file_ending = base_endings_r2_complete[base_endings_r1_complete.index(R1_file_ending)]
	fix_input_files_R2(base_endings_r2,".gz$",gzfiles,R2_file_ending)
elif len(gzfiles)==0 :
	base_endings_r1_complete = base_endings_r1
	base_endings_r2_complete = base_endings_r2
	R1_file_ending = fix_input_files_R1(base_endings_r1,"",onlyfiles,"TRUE")
	R2_file_ending = base_endings_r2_complete[base_endings_r1_complete.index(R1_file_ending)]
	fix_input_files_R2(base_endings_r2,"",onlyfiles,R2_file_ending)

##Convert any remaining categories
onlyfiles = [f for f in listdir(myfastqpath) if isfile(join(myfastqpath, f))]
gzfiles = list(filter(lambda x: re.search(r'.gz$', x), onlyfiles))
base_endings_r1_complete = base_endings_r1 + [s + ".gz" for s in base_endings_r1]
base_endings_r2_complete = base_endings_r2 + [s + ".gz" for s in base_endings_r2]
print(R1_file_ending)
print(R2_file_ending)

if len(gzfiles)>0 :
	# get files that aren't matching current file_ending
	odd_files = [i for i in [i for i in onlyfiles if not re.search(R1_file_ending, i)] if not re.search(R2_file_ending, i)]
	fastq_odd_1 = [i for i in odd_files if re.search(".fastq.gz", i)]
	fastq_odd_2 = [i for i in odd_files if re.search(".fq.gz", i)]
	if len(gzfiles) == len(odd_files) :
		print("Your dataset appears to be entirely single-end files.")
		fastq_odd = fastq_odd_1 + fastq_odd_2
		for myseq in range(len(fastq_odd)) :
			print(fastq_odd[myseq])
			print([i.replace(".fastq.gz","fq.gz") for i in fastq_odd][myseq])
			os.rename(join(myfastqpath,fastq_odd[myseq]),join(myfastqpath,[i.replace(".fastq.gz","fq.gz") for i in fastq_odd][myseq]))
	else :
		if len(fastq_odd_1)>0 :
			for myseq in range(len(fastq_odd_1)) :
				print(fastq_odd_1[myseq])
				print([i.replace(".fastq.gz",R1_file_ending+".gz") for i in fastq_odd_1][myseq])
				os.rename(join(myfastqpath,fastq_odd_1[myseq]),join(myfastqpath,[i.replace(".fastq.gz",R1_file_ending+".gz") for i in fastq_odd_1][myseq]))
		if len(fastq_odd_2)>0 :
			for myseq in range(len(fastq_odd_2)) :
				print(fastq_odd_2[myseq])
				print([i.replace(".fq.gz",R1_file_ending+".gz") for i in fastq_odd_2][myseq])
				os.rename(join(myfastqpath,fastq_odd_2[myseq]),join(myfastqpath,[i.replace(".fq.gz",R1_file_ending+".gz") for i in fastq_odd_2][myseq]))

elif len(gzfiles)==0 :
	odd_files = [i for i in [i for i in onlyfiles if not re.search(R1_file_ending, i)] if not re.search(R2_file_ending, i)]
	fastq_odd_1 = [i for i in odd_files if re.search(".fastq", i)]
	fastq_odd_2 = [i for i in odd_files if re.search(".fq", i)]
	if len(gzfiles) == len(odd_files) :
		print("Your dataset appears to be entirely single-end files.")
		fastq_odd = fastq_odd_1 + fastq_odd_2
		for myseq in range(len(fastq_odd)) :
			print(fastq_odd[myseq])
			print([i.replace(".fastq","fq") for i in fastq_odd][myseq])
			os.rename(join(myfastqpath,fastq_odd[myseq]),join(myfastqpath,[i.replace(".fastq","fq") for i in fastq_odd][myseq]))
	else :
		if len(fastq_odd_1)>0 :
			for myseq in range(len(fastq_odd_1)) :
				print(fastq_odd_1[myseq])
				print([i.replace(".fastq",R1_file_ending) for i in fastq_odd_1][myseq])
				os.rename(join(myfastqpath,fastq_odd_1[myseq]),join(myfastqpath,[i.replace(".fastq",R1_file_ending) for i in fastq_odd_1][myseq]))
		if len(fastq_odd_2)>0 :
			for myseq in range(len(fastq_odd_2)) :
				print(fastq_odd_2[myseq])
				print([i.replace(".fq",R1_file_ending) for i in fastq_odd_2][myseq])
				os.rename(join(myfastqpath,fastq_odd_2[myseq]),join(myfastqpath,[i.replace(".fq",R1_file_ending) for i in fastq_odd_2][myseq]))

onlyfiles = [f for f in listdir(myfastqpath) if isfile(join(myfastqpath, f))]
gzfiles = list(filter(lambda x: re.search(r'.gz$', x), onlyfiles))
if len(gzfiles)>0 :
	len_r1 = len([i for i in onlyfiles if re.search(re.escape(R1_file_ending) + ".gz$", i)])
else :
	len_r1 = len([i for i in onlyfiles if re.search(re.escape(R1_file_ending) + "$", i)])

if config["type"] == "single" :
	print("You have chosen to use single-end reads\nRead pairing not being checked...")

if config["type"] == "paired" :
	if len_r1*2 != len(onlyfiles) :
		myinput="One or more samples do not have a read pair!\nOr you are using a non-numeric read assignment!\nIf using paired-end samples, please ensure each sample has read 1 and read 2 files\nAborting..."
		raise NameError(myinput) ##exception to break workflow

sample_string = myfastqpath + "{sample}" + file_ending
print(sample_string)
SAMPLES, = glob_wildcards(sample_string)
print(SAMPLES)
print(R1_file_ending)
print(R2_file_ending)

def unique_function(x):
	return list(dict.fromkeys(x))

onlyfiles = [f for f in listdir(myfastqpath) if isfile(join(myfastqpath, f))]
gzfiles = list(filter(lambda x: re.search(r'.gz$', x), onlyfiles))
if len(gzfiles)>0 :
	if config["experiment"] == "cutrun" :
		my_cr_files = [i for i in gzfiles if re.search(re.escape(R1_file_ending) + ".gz$", i)] ##r1 files
		try :
			dedup_lengths = unique_function([len(next(SeqIO.parse(gzip.open(join(myfastqpath,i),"rt"), "fastq")).seq) for i in my_cr_files])
		except :
			raise NameError("One of your fastq files may be empty\nNow aborting...")

		if len(dedup_lengths)>1 :
			raise NameError("Based on sampling the first read of each R1 fastq file, your cut&run files have different read lengths!\nRecorded lengths:"+ ("elements in the list are "+', '.join(['%.f']*len(dedup_lengths))) % tuple(dedup_lengths) + " base pairs" + "\nAborting...")
		else :
			print("Congratulations, your cut&run fastq files appear to have uniform sequence lengths!\nProceeding with a read length of " + format(dedup_lengths[0]))
			config["read_length"] = dedup_lengths[0]
			config["read_length_max"] = dedup_lengths[0]-1

elif len(gzfiles) == 0 :
	if config["experiment"] == "cutrun" :
		my_cr_files = [i for i in onlyfiles if re.search(re.escape(R1_file_ending) + "$", i)] ##r1 files
		try :
			dedup_lengths = unique_function([len(next(SeqIO.parse(join(myfastqpath,i), "fastq")).seq) for i in my_cr_files])
		except :
			raise NameError("One of your fastq files may be empty\nNow aborting...")

		if len(dedup_lengths)>1 :
			raise NameError("Based on sampling the first read of each R1 fastq file, your cut&run files have different read lengths!\nRecorded lengths:"+ ("elements in the list are "+', '.join(['%.f']*len(dedup_lengths))) % tuple(dedup_lengths) + " base pairs" + "\nAborting...")
		else :
			print("Congratulations, your cut&run fastq files appear to have uniform sequence lengths!\nProceeding with a read length of " + format(dedup_lengths[0]))
			config["read_length"] = dedup_lengths[0]
			config["read_length_max"] = dedup_lengths[0]-1

##rules for resolving conflicting config parameters
if config["to_multiqc"] == "TRUE" : 
	for x in ["to_base_alignment","to_dedup_alignment","to_dedup_chr_alignment","to_dedup_tdf","to_dedup_bed","to_dedup_chr_bw","to_multiqc"] :
		config[x]="TRUE"

if config["to_dedup_chr_bw"] == "TRUE" : 
	for x in ["to_base_alignment","to_dedup_alignment","to_dedup_chr_alignment","to_dedup_bed"] :
		config[x]="TRUE"

if config["to_dedup_tdf"] == "TRUE" : ##check these
	for x in ["to_base_alignment","to_dedup_alignment","to_dedup_chr_alignment"] :
		config[x]="TRUE"

if config["to_dedup_bed"] == "TRUE" : ##check these in dag
	for x in ["to_base_alignment","to_dedup_alignment","to_dedup_chr_alignment"] :
		config[x]="TRUE"

if config["to_dedup_chr_alignment"] == "TRUE" : ##check these
	for x in ["to_base_alignment","to_dedup_alignment"] :
		config[x]="TRUE"

if config["to_dedup_alignment"] == "TRUE" :
	for x in ["to_base_alignment"] :
		config[x]="TRUE"

ALL_FASTQ = expand('output/trim_fastq/{sample}_R1_001_trimmed.fq.gz',sample=SAMPLES)
ALL_BAM = expand('output/bam/{sample}.bam',sample=SAMPLES)
ALL_DEDUP_BAM = expand('output/bam/{sample}.unique.sorted.rmdup.bam',sample=SAMPLES)
ALL_DEDUP_CHR_BAM = expand('output/bam/{sample}.unique.sorted.rmdup.chr.bam',sample=SAMPLES)
ALL_BW = expand('output/bw/{sample}.unique.sorted.rmdup.chr.bw',sample=SAMPLES)
ALL_BED = expand('output/bed/{sample}.unique.sorted.rmdup.chr.bed',sample=SAMPLES)
ALL_TDF = expand('output/tdf/{sample}.unique.sorted.rmdup.tdf',sample=SAMPLES)
COUNTS_MATRIX = "output/htseq_counts_matrix.txt"
MULTIQC_REPORT = "output/multiqc_report.html"

##Assign rule all inputs for config parameters
if config["experiment"] == "rnaseq":
	if config["keep_fastq"] == "TRUE" :
		if config["to_multiqc"] == "TRUE" : 
			rule all :
				input: 
					COUNTS_MATRIX, MULTIQC_REPORT, ALL_FASTQ
		elif config["to_base_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_BAM, ALL_FASTQ
		elif config["to_dedup_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_DEDUP_BAM, ALL_FASTQ
		elif config["to_dedup_chr_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_DEDUP_CHR_BAM, ALL_FASTQ

	elif config["keep_fastq"] == "FALSE" :
		if config["to_multiqc"] == "TRUE" : 
			rule all :
				input: 
					COUNTS_MATRIX, MULTIQC_REPORT
		elif config["to_base_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_BAM
		elif config["to_dedup_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_DEDUP_BAM
		elif config["to_dedup_chr_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_DEDUP_CHR_BAM

elif config["experiment"] == "chipseq" or config["experiment"] == "cutrun" :
	if config["keep_fastq"] == "TRUE" :
		if config["to_multiqc"] == "TRUE" : 
			rule all :
				input: 
					ALL_BED, ALL_TDF, ALL_BW, MULTIQC_REPORT, ALL_FASTQ
		elif config["to_base_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_BAM, ALL_FASTQ
		elif config["to_dedup_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_DEDUP_BAM, ALL_FASTQ
		elif config["to_dedup_chr_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_DEDUP_CHR_BAM, ALL_FASTQ
		elif config["to_dedup_tdf"] == "TRUE" : 
			rule all :
				input: 
					ALL_TDF, ALL_FASTQ
		elif config["to_dedup_bed"] == "TRUE" : 
			rule all :
				input: 
					ALL_BED, ALL_FASTQ
		elif config["to_dedup_chr_bw"] == "TRUE" : 
			rule all :
				input: 
					ALL_BW, ALL_FASTQ

	elif config["keep_fastq"] == "FALSE" :
		if config["to_multiqc"] == "TRUE" : 
			rule all :
				input: 
					ALL_BED, ALL_TDF, ALL_BW, MULTIQC_REPORT
		elif config["to_base_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_BAM
		elif config["to_dedup_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_DEDUP_BAM
		elif config["to_dedup_chr_alignment"] == "TRUE" :
			rule all :
				input: 
					ALL_DEDUP_CHR_BAM
		elif config["to_dedup_tdf"] == "TRUE" : 
			rule all :
				input: 
					ALL_TDF
		elif config["to_dedup_bed"] == "TRUE" : 
			rule all :
				input: 
					ALL_BED
		elif config["to_dedup_chr_bw"] == "TRUE" : 
			rule all :
				input: 
					ALL_BW

print("\"fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0] + "\"")
print("\"fastq/{sample}" + expand("{ending}",ending=R2_file_ending)[0] + "\"")

if config["experiment"] == "cutrun" and config["type"] == "paired" :
	rule trim_fastq_fastqc:
		input:
			pair1 = ("fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0]),
			pair2 = ("fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0])
		output:
			trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_val_1.fq.gz",
			trimmed_pair2 = "output/trim_fastq/{sample}_R2_001_val_2.fq.gz",
			fastqc_zipfile1 = "output/fastqc/{sample}_R1_001_fastqc.zip",
			fastqc_zipfile2 = "output/fastqc/{sample}_R2_001_fastqc.zip"
		log:
			"output/logs/{sample}.trim_adapters.log"
		run:
			shell("trim_galore {input.pair1} {input.pair2} --paired -o ./output/trim_fastq")
			shell("fastqc {input.pair1} {input.pair2} -o ./fastqc")
	rule split_length_long:
		input:
			trimg_pair1 = "output/trim_fastq/{sample}_R1_001_val_1.fq.gz",
			trimg_pair2 = "output/trim_fastq/{sample}_R2_001_val_2.fq.gz"
		output:
			cut_r1_p1 = temp("output/trim_fastq/{sample}_t1_R1.len75.fastq"),
			cut_r2_p1 = temp("output/trim_fastq/{sample}_t1_R2.len75.fastq")
		params:
			read_length = config["read_length"]
		log:
			"output/logs/{sample}.split_length_keeplong.log"
		run:
			shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r1_p1} {input.trimg_pair1}"),
			shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r2_p1} {input.trimg_pair2}")

	rule split_length_short:
		input:
			trimg_pair1 = "output/trim_fastq/{sample}_R1_001_val_1.fq.gz",
			trimg_pair2 = "output/trim_fastq/{sample}_R2_001_val_2.fq.gz"
		output:
			cut_r1_p2 = temp("output/trim_fastq/{sample}_t1_R1.lt75.fastq"),
			cut_r2_p2 = temp("output/trim_fastq/{sample}_t1_R2.lt75.fastq")
		params:
			read_length = config["read_length_max"]
		log:
			"output/logs/{sample}.split_length_keepshort.log"
		run:
			shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r1_p2} {input.trimg_pair1}"),
			shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r2_p2} {input.trimg_pair2}")
	rule trim_long:
		input:
			cut_r1_p1 = "output/trim_fastq/{sample}_t1_R1.len75.fastq",
			cut_r2_p1 = "output/trim_fastq/{sample}_t1_R2.len75.fastq"
		output:
			cut_r1_p3 = temp("output/trim_fastq/{sample}_t1_R1.len75_trim.fastq"),
			cut_r2_p3 = temp("output/trim_fastq/{sample}_t1_R2.len75_trim.fastq")
		log:
			"output/logs/{sample}.trim_long.log"
		run:
			shell("cutadapt -u -6 -o {output.cut_r1_p3} {input.cut_r1_p1}"),
			shell("cutadapt -u -6 -o {output.cut_r2_p3} {input.cut_r2_p1}")

	rule combine_split_lengths:
		input:
			cut_r1_p3 = "output/trim_fastq/{sample}_t1_R1.len75_trim.fastq",
			cut_r2_p3 = "output/trim_fastq/{sample}_t1_R2.len75_trim.fastq",
			cut_r1_p2 = "output/trim_fastq/{sample}_t1_R1.lt75.fastq",
			cut_r2_p2 = "output/trim_fastq/{sample}_t1_R2.lt75.fastq"
		output:
			cut_r1_p4 = temp("output/trim_fastq/{sample}_t2_R1.fastq"),
			cut_r2_p4 = temp("output/trim_fastq/{sample}_t2_R2.fastq")
		run:
			shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > {output.cut_r1_p4}"),
			shell("cat {input.cut_r2_p3} {input.cut_r2_p2} > {output.cut_r2_p4}")

	if config["keep_fastq"] == "FALSE" :
		rule sort_combined_lengths:
			input:
				cut_r1_p4 = "output/trim_fastq/{sample}_t2_R1.fastq",
				cut_r2_p4 = "output/trim_fastq/{sample}_t2_R2.fastq"
			output:
				cut_r1_p5 = temp("output/trim_fastq/{sample}_R1_001_trimmed.fq.gz"),
				cut_r2_p5 = temp("output/trim_fastq/{sample}_R2_001_trimmed.fq.gz")
			run:
				shell("cat {input.cut_r1_p4} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\\n\" > output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq"),
				shell("cat {input.cut_r2_p4} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\\n\" > output/trim_fastq/{wildcards.sample}_R2_001_trimmed.fq"),
				shell("gzip output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq"),
				shell("gzip output/trim_fastq/{wildcards.sample}_R2_001_trimmed.fq")

	if config["keep_fastq"] == "TRUE" :
		rule sort_combined_lengths:
			input:
				cut_r1_p4 = "output/trim_fastq/{sample}_t2_R1.fastq",
				cut_r2_p4 = "output/trim_fastq/{sample}_t2_R2.fastq"
			output:
				cut_r1_p5 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz",
				cut_r2_p5 = "output/trim_fastq/{sample}_R2_001_trimmed.fq.gz"
			run:
				shell("cat {input.cut_r1_p4} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\\n\" > output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq"),
				shell("cat {input.cut_r2_p4} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\\n\" > output/trim_fastq/{wildcards.sample}_R2_001_trimmed.fq"),
				shell("gzip output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq"),
				shell("gzip output/trim_fastq/{wildcards.sample}_R2_001_trimmed.fq")

	if config["to_base_alignment"] == "TRUE" and config["keep_base_alignment"] == "TRUE" :
		rule bowtie2:
			input:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz",
				trimmed_pair2 = "output/trim_fastq/{sample}_R2_001_trimmed.fq.gz"
			params:
				index = config["bowtie2_index"]
			output:
				bam = "output/bam/{sample}.bam",
				bambai = "output/bam/{sample}.bam.bai"
			threads:
				config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			run:
				shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} 2> {log} > processed/bam/{wildcards.sample}.sam"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm processed/bam/{wildcards.sample}.sam"),
				shell("samtools index {output.bam}")


	elif config["to_base_alignment"] == "TRUE" and config["keep_base_alignment"] == "FALSE" :
		rule bowtie2:
			input:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz",
				trimmed_pair2 = "output/trim_fastq/{sample}_R2_001_trimmed.fq.gz"
			params:
				index = config["bowtie2_index"]
			output:
				bam = temp("output/bam/{sample}.bam"),
				bambai = temp("output/bam/{sample}.bam.bai")
			threads:
				config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			run:
				shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} 2> {log} > processed/bam/{wildcards.sample}.sam"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm processed/bam/{wildcards.sample}.sam"),
				shell("samtools index {output.bam}")


	rule sort_alignment:
		input:
			"output/bam/{sample}.bam"
		output:
			bam = temp("output/bam/{sample}.sorted.bam"),
			bambai = temp("output/bam/{sample}.sorted.bam.bai")
		threads:
			config["threads_for_alignment"]
		log:
			"output/logs/{sample}.alignment.log"
		run:
			shell("samtools view -bh -f 3 -F 4 -F 8 {input} > output/bam/{wildcards.sample}_mapped.bam"),
			shell("samtools index output/bam/{wildcards.sample}_mapped.bam"),
			shell("samtools sort output/bam/{wildcards.sample}_mapped.bam > {output.bam}"),
			shell("samtools index {output.bam}")
			shell("rm output/bam/{wildcards.sample}_mapped.bam*")

	rule rmdup:
		input:
			"output/bam/{sample}.sorted.bam"
		output:
			bam = "output/bam/{sample}.unique.sorted.rmdup.bam"
		params:
			picardmetric = "output/logs/{sample}.markdups.metrics.txt"
		run:
			shell("picard MarkDuplicates INPUT={input} OUTPUT={output.bam} VALIDATION_STRINGENCY=SILENT METRICS_FILE={params.picardmetric}"),
			shell("samtools index {output.bam}")


if config["experiment"] == "cutrun" and config["type"] == "single" :
	rule trim_fastq_fastqc:
		input:
			pair1 = ("fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0])
		output:
			trimmed_pair1 = temp("output/trim_fastq/{sample}_R1_001_val_1.fq.gz"),
			fastqc_zipfile1 = "fastqc/{sample}_R1_001_fastqc.zip"
		log:
			"output/logs/{sample}.trim_adapters.log"
		run:
			shell("trim_galore {input.pair1} --paired -o ./output/trim_fastq")
			shell("fastqc {input.pair1} -o ./fastqc")

	rule split_length_long:
		input:
			trimg_pair1 = "output/trim_fastq/{sample}_R1_001_val_1.fq.gz"
		output:
			cut_r1_p1 = temp("output/trim_fastq/{sample}_t1_R1.len75.fastq")
		params:
			read_length = config["read_length"]
		log:
			"output/logs/{sample}.split_length_keeplong.log"
		run:
			shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r1_p1} {input.trimg_pair1}")
	
	rule split_length_short:
		input:
			trimg_pair1 = "output/trim_fastq/{sample}_R1_001_val_1.fq.gz"
		output:
			cut_r1_p2 = temp("output/trim_fastq/{sample}_t1_R1.lt75.fastq")
		params:
			read_length = config["read_length_max"]
		log:
			"output/logs/{sample}.split_length_keepshort.log"
		run:
			shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r1_p2} {input.trimg_pair1}")

	rule trim_long:
		input:
			cut_r1_p1 = "output/trim_fastq/{sample}_t1_R1.len75.fastq"
		output:
			cut_r1_p3 = temp("output/trim_fastq/{sample}_t1_R1.len75_trim.fastq")
		log:
			"output/logs/{sample}.trim_long.log"
		run:
			shell("cutadapt -u -6 -o {output.cut_r1_p3} {input.cut_r1_p1}")

	if config["keep_fastq"] == "FALSE" :
		rule combine_split_lengths:
			input:
				cut_r1_p3 = "output/trim_fastq/{sample}_t1_R1.len75_trim.fastq",
				cut_r1_p2 = "output/trim_fastq/{sample}_t1_R1.lt75.fastq"
			output:
				cut_r1_p4 = temp("output/trim_fastq/{sample}_R1_001_trimmed.fq.gz")
			run:
				shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > output/trim_fastq/{wildcards.sample}_R1_001_trimmed")
				shell("gzip output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq")

	elif config["keep_fastq"] == "TRUE" :
		rule combine_split_lengths:
			input:
				cut_r1_p3 = "output/trim_fastq/{sample}_t1_R1.len75_trim.fastq",
				cut_r1_p2 = "output/trim_fastq/{sample}_t1_R1.lt75.fastq"
			output:
				cut_r1_p4 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz"
			run:
				shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > output/trim_fastq/{wildcards.sample}_R1_001_trimmed")
				shell("gzip output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq")

	if config["to_base_alignment"] == "TRUE" and config["keep_base_alignment"] == "TRUE" :
		rule bowtie2:
			input:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz"
			params:
				index = config["bowtie2_index"]
			output:
				bam = "output/bam/{sample}.bam",
				bambai = "output/bam/{sample}.bam.bai"
			threads:
				config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			run:
				shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -U {input.trimmed_pair1} 2> {log} > output/bam/{wildcards.sample}.sam"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm output/bam/{wildcards.sample}.sam"),
				shell("samtools index {output.bam}")

	elif config["to_base_alignment"] == "TRUE" and config["keep_base_alignment"] == "FALSE" :
		rule bowtie2:
			input:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz"
			params:
				index = config["bowtie2_index"]
			output:
				bam = temp("output/bam/{sample}.bam"),
				bambai = temp("output/bam/{sample}.bam.bai")
			threads:
				config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			run:
				shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -U {input.trimmed_pair1} 2> {log} > output/bam/{wildcards.sample}.sam"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm output/bam/{wildcards.sample}.sam"),
				shell("samtools index {output.bam}")

	rule sort_alignment:
		input:
			"output/bam/{sample}.bam"
		output:
			bam = temp("output/bam/{sample}.sorted.bam"),
			bambai = temp("output/bam/{sample}.sorted.bam.bai")
		threads:
			config["threads_for_alignment"]
		log:
			"output/logs/{sample}.alignment.log"
		run:
			shell("samtools view -bh -f 3 -F 4 -F 8 {input} > output/bam/{wildcards.sample}_mapped.bam"),
			shell("samtools index output/bam/{wildcards.sample}_mapped.bam"),
			shell("samtools sort output/bam/{wildcards.sample}_mapped.bam > {output.bam}"),
			shell("samtools index {output.bam}")
			shell("rm output/bam/{wildcards.sample}_mapped.bam*")

	rule rmdup:
		input:
			"output/bam/{sample}.sorted.bam"
		output:
			bam = "output/bam/{sample}.unique.sorted.rmdup.bam"
		params:
			picardmetric = "output/logs/{sample}.markdups.metrics.txt"
		run:
			shell("picard MarkDuplicates INPUT={input} OUTPUT={output.bam} VALIDATION_STRINGENCY=SILENT METRICS_FILE={params.picardmetric}"),
			shell("samtools index {output.bam}")

if config["type"] == "single" and config["experiment"] != "cutrun" : # alignment
	if config["keep_fastq"] == "FALSE" :
		rule trim_fastq_fastqc:
			input:
				fastq = ("fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0])
			output:
				trimmed_fastq = temp("output/trim_fastq/{sample}_R1_001_trimmed.fq.gz"),
				fastqc_zipfile = "output/fastqc/{sample}_R1_001_fastqc.zip"
			log:
				"output/logs/{sample}.trim_adapters.log"
			run:
				shell("trim_galore {input.fastq} -o ./output/trim_fastq")
				shell("fastqc {input.fastq} -o ./output/fastqc")		
	if config["keep_fastq"] == "TRUE" :
		rule trim_fastq_fastqc:
			input:
				fastq = ("fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0])
			output:
				trimmed_fastq = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz",
				fastqc_zipfile = "output/fastqc/{sample}_R1_001_fastqc.zip"
			log:
				"output/logs/{sample}.trim_adapters.log"
			run:
				shell("trim_galore {input.fastq} -o ./output/trim_fastq")
				shell("fastqc {input.fastq} -o ./output/fastqc")
	if config["to_base_alignment"] == "TRUE" and config["keep_base_alignment"] == "TRUE" :
		rule fastq_to_sam:
			input:
				trimmed_fastq = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz"
			params:
				index = config["index"]
			output:
				bam = "output/bam/{sample}.bam",
				bambai = "output/bam/{sample}.bam.bai"
			threads: config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			run:
				shell("hisat2 -p {threads} -x {params.index} -U {input.trimmed_fastq} -S {output.sam} 2> {log}"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm output/bam/{wildcards.sample}.sam"),
				shell("samtools index {output.bam}")

	if config["to_base_alignment"] == "TRUE" and config["keep_base_alignment"] == "FALSE" :
		rule fastq_to_sam:
			input:
				trimmed_fastq = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz"
			params:
				index = config["index"]
			output:
				bam = temp("output/bam/{sample}.bam")
			threads: config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			run:
				shell("hisat2 -p {threads} -x {params.index} -U {input.trimmed_fastq} -S {output.sam} 2> {log}")
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm output/bam/{wildcards.sample}.sam")


elif config["type"] == "paired" and config["experiment"] != "cutrun" :
	if config["keep_fastq"] == "FALSE" :
		rule trim_fastq_fastqc:
			input:
		                pair1 = ("fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0]),
		                pair2 = ("fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0])
			output:
				trimmed_pair1 = temp("output/trim_fastq/{sample}_R1_001_val_1.fq.gz"),
				trimmed_pair2 = temp("output/trim_fastq/{sample}_R2_001_val_2.fq.gz"),
				fastqc_zipfile1 = "fastqc/{sample}_R1_001_fastqc.zip",
				fastqc_zipfile2 = "fastqc/{sample}_R2_001_fastqc.zip"
			log:
				"output/logs/{sample}.trim_adapters.log"
			run:
				shell("trim_galore {input.pair1} {input.pair2} --paired -o ./output/trim_fastq")
				shell("fastqc {input.pair1} {input.pair2} -o ./fastqc")
	if config["keep_fastq"] == "TRUE" :
		rule trim_fastq_fastqc:
			input:
		                pair1 = ("fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0]),
		                pair2 = ("fastq/{sample}" + expand("{ending}",ending=R1_file_ending)[0])
			output:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_val_1.fq.gz",
				trimmed_pair2 = "output/trim_fastq/{sample}_R2_001_val_2.fq.gz",
				fastqc_zipfile1 = "fastqc/{sample}_R1_001_fastqc.zip",
				fastqc_zipfile2 = "fastqc/{sample}_R2_001_fastqc.zip"
			log:
				"output/logs/{sample}.trim_adapters.log"
			run:
				shell("trim_galore {input.pair1} {input.pair2} --paired -o ./output/trim_fastq")
				shell("fastqc {input.pair1} {input.pair2} -o ./fastqc")

	if config["to_base_alignment"] == "TRUE" and config["keep_base_alignment"] == "TRUE" :
		rule fastq_to_sam:
			input:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_val_1.fq.gz",
				trimmed_pair2 = "output/trim_fastq/{sample}_R2_001_val_2.fq.gz"
			params:
				index = config["index"]
			output:
				bam = "output/bam/{sample}.bam",
				bambai = "output/bam/{sample}.bam.bai"
			threads: config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			shell:
				shell("hisat2 -p {threads} -x {params.index} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} -S {output.sam} 2> {log}"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm output/bam/{wildcards.sample}.sam"),
				shell("samtools index {output.bam}")

	if config["to_base_alignment"] == "TRUE" and config["keep_base_alignment"] == "FALSE" :
		rule fastq_to_bam:
			input:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_val_1.fq.gz",
				trimmed_pair2 = "output/trim_fastq/{sample}_R2_001_val_2.fq.gz"
			params:
				index = config["index"]
			output:
				bam = temp("output/bam/{sample}.bam"),
				bambai = temp("output/bam/{sample}.bam.bai")
			threads: config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			run:
				shell("hisat2 -p {threads} -x {params.index} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} -S {output.sam} 2> {log}"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm output/bam/{wildcards.sample}.sam")

if config["experiment"] == "chipseq" :
	rule bam_to_unique_mapped:	##Remove the multimapped reads and sort
		input:
			"output/bam/{sample}.bam"
		output:
			bam = temp("output/bam/{sample}.sorted.bam")
		run:
			shell("samtools view -bh -f 3 -F 4 -F 8 -F 256 {input} > output/bam/{wildcards.sample}_filtered.bam") ##mapped, unique
			shell("samtools sort -O BAM -o {output.bam} output/bam/{wildcards.sample}_filtered.bam")
			shell("rm output/bam/{wildcards.sample}_filtered.bam")

if config["experiment"] == "rnaseq":
	rule sortedbam_to_rmdup_rnaseq:
		input:
			sorted_bam = "output/bam/{sample}.bam"
		output:
			dup_removed = "output/bam/{sample}.sorted.rmdup.bam"
		log:
			"output/logs/{sample}.rmdup.log"
		run:
			shell("samtools rmdup {input.sorted_bam} {output.dup_removed} 2> {log}")
			shell("rm -f {input.sorted_bam}")

# chipseq 
if config["experiment"] == "chipseq":
	rule sortedbam_to_rmdup:
		input:
			sorted_bam = "output/bam/{sample}.sorted.bam"
		output:
			dup_removed = "output/bam/{sample}.unique.sorted.rmdup.bam"
		log:
			"output/logs/{sample}.rmdup.log"
		run:
			shell("samtools rmdup {input.sorted_bam} {output.dup_removed} 2> {log}")
			shell("rm -f {input.sorted_bam}")

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
		chrbam = "output/bam/{sample}.unique.sorted.rmdup.chr.bam"
	log:
		"output/logs/{sample}.chrbam.log"
	shell:
		'samtools view -H {input.dup_removed} '
		'| sed -e "s/SN:\([0-9XY]\)/SN:chr\\1/" -e "s/SN:MT/SN:chrM/" '
		'| samtools reheader - {input.dup_removed} > {output.chrbam}'

rule chrbam_to_bw:
	input:
		chrbam = "output/bam/{sample}.unique.sorted.rmdup.chr.bam"
	output:
		bw_file = "output/bw/{sample}.unique.sorted.rmdup.chr.bw"
	log:
		"output/logs/{sample}.bw.log"
	run:
		shell("samtools index {input.chrbam}")
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
if config["count_scheme"] == "fraction":
	rule sortedbam_to_counts:
		input:
			sorted_bam = "output/bam/{sample}.sorted.bam" if config["type"] == "single"
							else "output/bam/{sample}.sorted.rmdup.bam",
			gtf = config["gtf"]
		output:
			counts = "output/counts/{sample}.counts.txt"
		log:
			"output/logs/{sample}.htseq_counts.log"
		shell:
			"featureCounts -p -O --fraction -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}" if config["type"] == "paired" 
			else "featureCounts -O --fraction -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}"

elif config["count_scheme"] == "count_all":
	rule sortedbam_to_counts:
		input:
			sorted_bam = "output/bam/{sample}.sorted.bam" if config["type"] == "single"
							else "output/bam/{sample}.sorted.rmdup.bam",
			gtf = config["gtf"]
		output:
			counts = "output/counts/{sample}.counts.txt"
		log:
			"output/logs/{sample}.htseq_counts.log"
		shell:
			"featureCounts -p -O -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}" if config["type"] == "paired" 
			else "featureCounts -O -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}"

elif config["count_scheme"] == "count_uniq":
	rule sortedbam_to_counts:
		input:
			sorted_bam = "output/bam/{sample}.sorted.bam" if config["type"] == "single"
							else "output/bam/{sample}.sorted.rmdup.bam",
			gtf = config["gtf"]
		output:
			counts = "output/counts/{sample}.counts.txt"
		log:
			"output/logs/{sample}.htseq_counts.log"
		shell:
			"featureCounts -p -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}" if config["type"] == "paired" 
			else "featureCounts -t gene -a {input.gtf} -o {output.counts} "
			"{input.sorted_bam} 2> {log}"
			

rule counts_matrix:
	input:
		counts = expand("output/counts/{sample}.counts.txt", \
						sample=SAMPLES)
	output:
		matrix = "output/htseq_counts_matrix.txt"
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
			matrix = "output/htseq_counts_matrix.txt"
		output:
			multiqc_report = "multiqc_report.html"
		params:
			multiqc_config = config["multiqc_yaml"]
		shell:
			"multiqc . -f --config {params.multiqc_config}"

elif config["experiment"] == "chipseq" or config["experiment"] == "cutrun" :
	rule run_multiqc:
		input:
			tdf = expand("output/tdf/{sample}.unique.sorted.rmdup.tdf", \
						sample = SAMPLES)
		output:
			multiqc_report = "multiqc_report.html"
		params:
			multiqc_config = config["multiqc_yaml"]
		shell:
			"multiqc . -f --config {params.multiqc_config}"



