configfile: "config.yaml"

myfastqpath="./fastq/"
import os
from os import listdir
from os.path import isfile, join
import re
from Bio import SeqIO
import gzip
##List the files present in the fastq folder.
onlyfiles = [f for f in listdir(myfastqpath) if isfile(join(myfastqpath, f))]
l = [list(filter(lambda x: re.search(i, x), onlyfiles)) for i in [".fastq$",".fq$",".fastq.gz$",".fq.gz$"]]
onlyfiles = [item for sublist in l for item in sublist]
gzfiles = list(filter(lambda x: re.search(r'.gz$', x), onlyfiles))

if len(gzfiles)==0 and len(onlyfiles)==0 :
	raise NameError("You do not seem to have any fastq files present to process. Now exiting...") ##exception to break workflow


if len(gzfiles)>0 :
	if len(gzfiles)!=len(onlyfiles) : 	##Error message:
		myinput="You have a mixture of gzipped files and non-gzipped files\nOnly {} of total {} files are gzipped!"
		raise NameError(print(myinput.format(len(gzfiles),len(onlyfiles)))) ##exception to break workflow

##Create the pattern of file endings
def create_endings(x) :
	return(["_R"+format(x)+"_001.fastq", "_R"+format(x)+"_001.fq", "_R"+format(x)+".fastq", "_R"+format(x)+".fq", "_"+format(x)+".fastq","_"+format(x)+".fq", ".R"+format(x)+"_001.fastq", ".R"+format(x)+"_001.fq", ".R"+format(x)+".fastq", ".R"+format(x)+".fq", "."+format(x)+".fastq", "."+format(x)+".fq","_r"+format(x)+"_001.fastq", "_r"+format(x)+"_001.fq", "_r"+format(x)+".fastq", "_r"+format(x)+".fq", ".r"+format(x)+"_001.fastq", ".r"+format(x)+"_001.fq", ".r"+format(x)+".fastq", ".r"+format(x)+".fq"])

base_endings_r1 = create_endings(1)
base_endings_r2 = create_endings(2)

def fix_input_files(file_suffix,input_fileset) :
	if file_suffix == ".gz" :
		ending_dictionary = dict(zip(base_endings_r1,base_endings_r2))
	else :
		ending_dictionary = dict(zip(base_endings_r1,base_endings_r2)) ##Define the R1 and R1 suffix pairs for reference
	mylist = list() ##first traverse the R1 base endings to find the common ending	
	for x in base_endings_r1 :
		my_regex = re.escape(x) + file_suffix
		if(len([i for i in input_fileset if re.search(my_regex, i)])>0) :
			mylist.append((x,len([i for i in input_fileset if re.search(my_regex, i)])))
	mylist_dict = {key: value for (key, value) in mylist}
	inverse = [(value, key) for key, value in mylist_dict.items()]
	myR1_suffix = max(inverse)[1]
	myR2_suffix = ending_dictionary[myR1_suffix]
	for i in ("R1","R2") :
		if i=="R1" :
			suffixes_to_change = [x for x in list(mylist_dict.keys()) if x!=myR1_suffix]
			if len(suffixes_to_change)>0 :
				for x in suffixes_to_change :
					my_regex = re.escape(x) + file_suffix
					oldnames = [i for i in input_fileset if re.search(my_regex, i)]
					print(oldnames) ##
					for myseq in range(len(oldnames)) :
						print(oldnames[myseq]) ##
						print([i.replace(x,myR1_suffix) for i in oldnames][myseq]) ##
						os.rename(join(myfastqpath,oldnames[myseq]),join(myfastqpath,[i.replace(x,myR1_suffix) for i in oldnames][myseq]))
		if i=="R2" :
			mylist = list()
			for x in base_endings_r2 :
				my_regex = re.escape(x) + file_suffix
				if(len([i for i in input_fileset if re.search(my_regex, i)])>0) :
					mylist.append((x,len([i for i in input_fileset if re.search(my_regex, i)])))
			suffixes_to_change = [x for x,y in mylist if x!=myR2_suffix]
			if len(suffixes_to_change)>0 :
				for x in suffixes_to_change :
					my_regex = re.escape(x) + file_suffix
					oldnames = [i for i in input_fileset if re.search(my_regex, i)]
					print(oldnames) ##
					for myseq in range(len(oldnames)) :
						print(oldnames[myseq]) ##
						print([i.replace(x,myR2_suffix) for i in oldnames][myseq]) ##
						os.rename(join(myfastqpath,oldnames[myseq]),join(myfastqpath,[i.replace(x,myR2_suffix) for i in oldnames][myseq]))
	return([myR1_suffix,myR2_suffix])


if len(gzfiles)>0 :
	good_endings = fix_input_files(".gz",gzfiles)
else :
	good_endings = fix_input_files("",onlyfiles)

onlyfiles = [f for f in listdir(myfastqpath) if isfile(join(myfastqpath, f))]
l = [list(filter(lambda x: re.search(i, x), onlyfiles)) for i in [".fastq$",".fq$",".fastq.gz$",".fq.gz$"]]
onlyfiles = [item for sublist in l for item in sublist]
gzfiles = list(filter(lambda x: re.search(r'.gz$', x), onlyfiles))

def convert_single(suffix,input_fileset) :
	R1_file_ending = good_endings[0]
	R2_file_ending = good_endings[1]
	odd_files = [i for i in [i for i in input_fileset if not re.search(R1_file_ending, i)] if not re.search(R2_file_ending, i)]
	fastq_odd_1 = [i for i in odd_files if re.search(".fastq"+suffix+"$", i)]
	fastq_odd_2 = [i for i in odd_files if re.search(".fq"+suffix+"$", i)]
	if len(input_fileset) == len(odd_files) :
		print("Your dataset appears to be entirely single-end files.")
	if len(odd_files) > 0 :
		print("Now unifying "+format(len(odd_files))+" single-end files to \""+R1_file_ending+suffix+"\" ending")
		if len(fastq_odd_1)>0 :
			for myseq in range(len(fastq_odd_1)) :
				print(fastq_odd_1[myseq])
				print([i.replace(".fastq"+suffix,R1_file_ending+suffix) for i in fastq_odd_1][myseq])
				os.rename(join(myfastqpath,fastq_odd_1[myseq]),join(myfastqpath,[i.replace(".fastq"+suffix,R1_file_ending+suffix) for i in fastq_odd_1][myseq])) # rename files to correct ending
		if len(fastq_odd_2)>0 :
			for myseq in range(len(fastq_odd_2)) :
				print(fastq_odd_2[myseq])
				print([i.replace(".fq"+suffix,R1_file_ending+suffix) for i in fastq_odd_2][myseq])
				os.rename(join(myfastqpath,fastq_odd_2[myseq]),join(myfastqpath,[i.replace(".fq"+suffix,R1_file_ending+suffix) for i in fastq_odd_2][myseq])) # rename files to correct ending			

if len(gzfiles)>0 :
	convert_single(".gz",gzfiles)
else :
	convert_single("",onlyfiles)

onlyfiles = [f for f in listdir(myfastqpath) if isfile(join(myfastqpath, f))]
l = [list(filter(lambda x: re.search(i, x), onlyfiles)) for i in [".fastq$",".fq$",".fastq.gz$",".fq.gz$"]]
onlyfiles = [item for sublist in l for item in sublist]
gzfiles = list(filter(lambda x: re.search(r'.gz$', x), onlyfiles))

##Now check the file pairing
if config["type"] == "single" :
	print("You have chosen to use single-end reads\nRead pairing not being checked...")
elif config["type"] == "paired" :
	if len(gzfiles)>0 :
		len_r1 = len([i for i in onlyfiles if re.search(re.escape(good_endings[0]) + ".gz$", i)])
	else :
		len_r1 = len([i for i in onlyfiles if re.search(re.escape(good_endings[0]) + "$", i)])
	if len_r1*2 != len(onlyfiles) :
		myinput="One or more samples do not have a read pair!\nIf using paired-end samples, please ensure each sample has read 1 and read 2 files\nAborting..."
		raise NameError(myinput) ##exception to break workflow
else :
	myinput="You have specified read type: "+config["type"]+"\nPlease specify either \"paired\" or \"single\" in the config.yaml file, then rerun the pipeline."
	raise NameError(myinput)

R1_file_ending = good_endings[0]
R2_file_ending = good_endings[1]
if len(gzfiles)>0 :
	sample_string = myfastqpath + "{sample}" + R1_file_ending+".gz"
	suffix=".gz"
else :
	sample_string = myfastqpath + "{sample}" + R1_file_ending
	suffix=""

SAMPLES, = glob_wildcards(sample_string)

def unique_function(x) :
	return list(dict.fromkeys(x))
def check_readlength(suffix,input_fileset) :
	my_cr_files = [i for i in gzfiles if re.search(re.escape(R1_file_ending) + suffix + "$", i)] ##r1 files
	try :
		dedup_lengths = unique_function([len(next(SeqIO.parse(gzip.open(join(myfastqpath,i),"rt"), "fastq")).seq) for i in my_cr_files])
	except :
		raise NameError("One of your fastq files may be empty\nNow aborting...")
	if len(dedup_lengths)>1 :
		raise NameError("Based on sampling the first read of each R1 fastq file, your cut&run files have different read lengths!\nRecorded lengths:"+ ("elements in the list are "+', '.join(['%.f']*len(dedup_lengths))) % tuple(dedup_lengths) + " base pairs" + "\nAborting...")
	else :
		print("Congratulations, your cut&run fastq files appear to have uniform sequence lengths!\nProceeding with a read length of " + format(dedup_lengths[0]))
		return(dedup_lengths[0],dedup_lengths[0]-1)


if len(gzfiles) > 0 and config["experiment"] == "cutrun" :
	read_len_check = check_readlength(".gz",gzfiles)
	config["read_length"] = read_len_check[0]
	config["read_length_max"] = read_len_check[0]-1

elif len(gzfiles) == 0 and config["experiment"] == "cutrun" :
	check_readlength("",onlyfiles)
	config["read_length"] = read_len_check[0]
	config["read_length_max"] = read_len_check[0]-1

##Create function for creating rule sets
def choose_rule_all(config) :
	myout = []
	if config["to_multiqc"] == "TRUE" :
		myout.append("output/multiqc_report.html")
	if config["to_bw"] == "TRUE" and config["experiment"] != "rnaseq":
		myout.append(expand('output/bw/{sample}.unique.sorted.rmdup.chr.bw',sample=SAMPLES))
	if config["to_bed"] == "TRUE" and config["experiment"] != "rnaseq" :
		myout.append(expand('output/bed/{sample}.unique.sorted.rmdup.chr.bed',sample=SAMPLES))
	if config["to_tdf"] == "TRUE" and config["experiment"] != "rnaseq" :
		myout.append(expand('output/tdf/{sample}.unique.sorted.rmdup.tdf',sample=SAMPLES))
	if config["keep_unfiltered_bam"] == "TRUE" :
		myout.append(expand('output/bam/{sample}.bam',sample=SAMPLES))
	if config["keep_fastq"] == "TRUE" :
		myout.append(expand('output/trim_fastq/{sample}_R1_001_trimmed.fq.gz',sample=SAMPLES))
	if config["experiment"] == "rnaseq" :
		myout.append("output/htseq_counts_matrix.txt")
		myout.append(expand('output/fastqc/{sample}_R1_001_fastqc.zip',sample=SAMPLES))
	return(myout)

rule all :
	input:  
		choose_rule_all(config)

print(choose_rule_all(config))

def create_inputs(config) :
	return([("fastq/{sample}" + expand("{ending}{suffix}",ending=R1_file_ending,suffix=suffix)[0]+""),("fastq/{sample}" + expand("{ending}{suffix}",ending=R2_file_ending,suffix=suffix)[0]+"")])

if config["experiment"] == "cutrun" :
	if config["type"] == "paired" :
		rule trim_fastq_fastqc :
			input :
				pair1 = create_inputs(config)[0],
				pair2 = create_inputs(config)[1]
			output:
				trimmed_pair1 = temp("output/trim_fastq/{sample}_R1_001_val_1.fq.gz"),
				trimmed_pair2 = temp("output/trim_fastq/{sample}_R2_001_val_2.fq.gz"),
				fastqc_zipfile1 = "output/fastqc/{sample}_R1_001_fastqc.zip",
				fastqc_zipfile2 = "output/fastqc/{sample}_R2_001_fastqc.zip" 
			log:
				"output/logs/{sample}.trim_adapters.log"
			run:
				shell("trim_galore {input.pair1} {input.pair2} --paired -o ./output/trim_fastq")
				shell("fastqc {input.pair1} {input.pair2} -o ./output/fastqc")

	if config["type"] == "single" :
		rule trim_fastq_fastqc :
			input :
				pair1 = create_inputs(config)[0]
			output:
				trimmed_pair1 = temp("output/trim_fastq/{sample}_R1_001_val_1.fq.gz"),
				trimmed_pair2 = temp("output/trim_fastq/{sample}_R2_001_val_2.fq.gz"),
				fastqc_zipfile1 = "output/fastqc/{sample}_R1_001_fastqc.zip",
				fastqc_zipfile2 = "output/fastqc/{sample}_R2_001_fastqc.zip" 
			log:
				"output/logs/{sample}.trim_adapters.log"
			run:
				shell("trim_galore {input.pair1} -o ./output/trim_fastq --basename Temp_R1_{wildcards.sample}")
				shell("mv output/trim_fastq/Temp_R1_{wildcards.sample}_trimmed.fq.gz output/trim_fastq/{wildcards.sample}_R1_001_val_1.fq.gz")
				shell("fastqc {input.pair1} -o ./output/fastqc")
				shell("touch {output.trimmed_pair2}")
				shell("touch {output.fastqc_zipfile2}")

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
			if config["type"] == "paired" :
				shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r1_p1} {input.trimg_pair1}"),
				shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r2_p1} {input.trimg_pair2}")
			else :
				shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r1_p1} {input.trimg_pair1}")
				shell("touch {output.cut_r2_p1}")
				shell("rm output/fastqc/{sample}_R2_001_fastqc.zip") 

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
			if config["type"] == "paired" :
				shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r1_p2} {input.trimg_pair1}"),
				shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r2_p2} {input.trimg_pair2}")
			else :
				shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r1_p2} {input}")
				shell("touch {output.cut_r2_p2}")

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
			if config["type"] == "paired" :
				shell("cutadapt -u -6 -o {output.cut_r1_p3} {input.cut_r1_p1}"),
				shell("cutadapt -u -6 -o {output.cut_r2_p3} {input.cut_r2_p1}")
			else :
				shell("cutadapt -u -6 -o {output.cut_r1_p3} {input}")
				shell("touch {output.cut_r2_p3}")


	rule combine_split_lengths:
		input:
			cut_r1_p3 = "output/trim_fastq/{sample}_t1_R1.len75_trim.fastq",
			cut_r2_p3 = "output/trim_fastq/{sample}_t1_R2.len75_trim.fastq",
			cut_r1_p2 = "output/trim_fastq/{sample}_t1_R1.lt75.fastq",
			cut_r2_p2 = "output/trim_fastq/{sample}_t1_R2.lt75.fastq"
		output:
			cut_r1_p4 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz",
			cut_r2_p4 = "output/trim_fastq/{sample}_R2_001_trimmed.fq.gz"		
		run:
			if config["type"] == "paired" :
				shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > output/trim_fastq/{wildcards.sample}_t2_R1.fastq"),
				shell("cat {input.cut_r2_p3} {input.cut_r2_p2} > output/trim_fastq/{wildcards.sample}_t2_R2.fastq"),
				shell("cat output/trim_fastq/{wildcards.sample}_t2_R1.fastq | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\\n\" > output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq"),
				shell("cat output/trim_fastq/{wildcards.sample}_t2_R2.fastq | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\\n\" > output/trim_fastq/{wildcards.sample}_R2_001_trimmed.fq"),
				shell("gzip output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq"),
				shell("gzip output/trim_fastq/{wildcards.sample}_R2_001_trimmed.fq")
			else :
				shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > output/trim_fastq/{wildcards.sample}_R1_001_trimmed")
				shell("gzip output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq")
				shell("touch {output.cut_r2_p4}")

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
			if config["type"] == "paired" :
				shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} 2> {log} > output/bam/{wildcards.sample}.sam"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm output/bam/{wildcards.sample}.sam"),
				shell("samtools index {output.bam}")
			else : 
				shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -U {input.trimmed_pair1} 2> {log} > output/bam/{wildcards.sample}.sam"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm output/bam/{wildcards.sample}.sam"),
				shell("samtools index {output.bam}")
				shell("rm {input.trimmed_pair2}")
			if config["keep_fastq"] == "FALSE" and config["type"] == "paired" :
				shell("rm {input.trimmed_pair1} {input.trimmed_pair2}")
			elif config["keep_fastq"] == "FALSE" and config["type"] == "single" :
				shell("rm {input.trimmed_pair1}")

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
			if config["keep_unfiltered_bam"] == "FALSE" :
				shell("rm output/bam/{wildcards.sample}.bam")
				shell("rm output/bam/{wildcards.sample}.bam.bai")

	rule rmdup:
		input:
			"output/bam/{sample}.sorted.bam"
		output:
			bam = "output/bam/{sample}.unique.sorted.rmdup.bam",
			bambai = "output/bam/{sample}.unique.sorted.rmdup.bam.bai"
		params:
			picardmetric = "output/logs/{sample}.markdups.metrics.txt"
		run:
			shell("picard MarkDuplicates INPUT={input} OUTPUT={output.bam} VALIDATION_STRINGENCY=SILENT METRICS_FILE={params.picardmetric}"),
			shell("samtools index {output.bam}")



def create_inputs(config) :
	return([("fastq/{sample}" + expand("{ending}{suffix}",ending=R1_file_ending,suffix=suffix)[0]+""),("fastq/{sample}" + expand("{ending}{suffix}",ending=R2_file_ending,suffix=suffix)[0]+"")])

if config["experiment"] == "rnaseq" or config["experiment"] == "chipseq" : 
	if config["type"] == "paired" :
		rule trim_fastq_fastqc:
			input :
				pair1 = create_inputs(config)[0],
				pair2 = create_inputs(config)[1]
			output:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz",
				trimmed_pair2 = "output/trim_fastq/{sample}_R2_001_trimmed.fq.gz",
				fastqc_zipfile1 = "output/fastqc/{sample}_R1_001_fastqc.zip",
				fastqc_zipfile2 = "output/fastqc/{sample}_R2_001_fastqc.zip"			
			log:
				"output/logs/{sample}.trim_adapters.log"
			run:
				shell("trim_galore {input.pair1} {input.pair2} --paired -o ./output/trim_fastq"),
				shell("fastqc {input.pair1} {input.pair2} -o ./output/fastqc"),
				shell("mv output/trim_fastq/{wildcards.sample}_R1_001_val_1.fq.gz output/trim_fastq/{wildcards.sample}_R1_001_trimmed.fq.gz"),
				shell("mv output/trim_fastq/{wildcards.sample}_R2_001_val_2.fq.gz output/trim_fastq/{wildcards.sample}_R2_001_trimmed.fq.gz")

		rule fastq_to_bam:
			input:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz",
				trimmed_pair2 = "output/trim_fastq/{sample}_R2_001_trimmed.fq.gz"
			params:
				index = config["index"]
			output:
				bam = "output/bam/{sample}.bam",
				bambai = "output/bam/{sample}.bam.bai"
			threads: config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			run:
				shell("hisat2 -p {threads} -x {params.index} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} -S output/bam/{wildcards.sample}.sam 2> {log}"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("samtools index {output.bam}")
				shell("rm output/bam/{wildcards.sample}.sam")
				if config["keep_fastq"] == "FALSE" :
					shell("rm {input.trimmed_pair1} {input.trimmed_pair2}")

	elif config["type"] == "single" :
		print(create_inputs(config)[0])
		rule trim_fastq_fastqc:
			input :
				pair1 = create_inputs(config)[0]
			output:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz",
				fastqc_zipfile1 = "output/fastqc/{sample}_R1_001_fastqc.zip"		
			log:
				"output/logs/{sample}.trim_adapters.log"
			run:
				shell("trim_galore {input.pair1} -o ./output/trim_fastq"),
				shell("fastqc {input.pair1} -o ./output/fastqc")

		rule fastq_to_bam:
			input:
				trimmed_pair1 = "output/trim_fastq/{sample}_R1_001_trimmed.fq.gz"
			params:
				index = config["index"]
			output:
				bam = "output/bam/{sample}.bam",
				bambai = "output/bam/{sample}.bam.bai"
			threads: config["threads_for_alignment"]
			log:
				"output/logs/{sample}.alignment.log"
			run:
				shell("hisat2 -p {threads} -x {params.index} -U {input.trimmed_pair1} -S output/bam/{wildcards.sample}.sam 2> {log}"),
				shell("samtools sort output/bam/{wildcards.sample}.sam | samtools view -bS - > {output.bam}"),
				shell("rm output/bam/{wildcards.sample}.sam"),
				shell("samtools index {output.bam}")
				if config["keep_fastq"] == "FALSE" :
					shell("rm {input.trimmed_pair1}")


	if config["experiment"] == "chipseq" :
		rule bam_to_unique_mapped:	##Remove the multimapped reads and sort
			input:
				"output/bam/{sample}.bam"
			output:
				bam = temp("output/bam/{sample}.sorted.bam")
			run:
				if config["type"] == "paired" :
					shell("samtools view -bh -f 3 -F 4 -F 8 -F 256 {input} > output/bam/{wildcards.sample}_filtered.bam") ##mapped with pair, unique
					shell("samtools sort -O BAM -o {output.bam} output/bam/{wildcards.sample}_filtered.bam")
					shell("rm output/bam/{wildcards.sample}_filtered.bam")
				else : 
					shell("samtools view -bh -F 4 -F 256 {input} > output/bam/{wildcards.sample}_filtered.bam") ##mapped with pair, unique locations
					shell("samtools sort -O BAM -o {output.bam} output/bam/{wildcards.sample}_filtered.bam")
					shell("rm output/bam/{wildcards.sample}_filtered.bam")
				if config["keep_unfiltered_bam"] == "FALSE" :
					shell("rm output/bam/{sample}.bam output/bam/{sample}.bam.bai")

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
				"output/bam/{sample}.bam"
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
				"output/bam/{sample}.bam"
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
