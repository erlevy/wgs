# UC San Diego Medical School - Division of Biomedical Informatics
#
# Author: Hannah Burkhardt
# Date: 7/18/2013
#
# Filename: kd_pipeline.py
# Purpose: To generate a vcf file with SNPs from fastq raw reads.
#
#


#####################################################
##--import statements, usage statement, read args--##
#####################################################

import sys
import subprocess
from datetime import datetime
from time import time

# print usage statement
if len(sys.argv) < 8:
	print "--- usage: python kd_pipeline.py <GATK_path> <Picard_path> <snpEff_path> <ref.fa> <reads1.fastq> <reads2.fastq> <known_sites.vcf> <create new folder? 0/1>"
	print "version: 3.0"
	print "Improved in this version:\n- aln/sampe.\n- output files in own subdirectory\n- Validation Stringency is lenient\n- Readgroups added with AddOrReplaceReadGroups rather than bwa\n- automatic logfile in addition to stdout log\n- Mulitple performance options, tailored to each tool\n- file sizes printed"
	sys.exit(1)

start_time = datetime.now()

gatk = sys.argv[1]
picard = sys.argv[2]
snpeff = sys.argv[3]
ref = sys.argv[4]
reads1 = sys.argv[5]
reads2 = sys.argv[6]
knownsites = sys.argv[7]
timestamp =  str(start_time.date()) + "__" + str(start_time.hour) + "-" + str(start_time.minute)
comment = sys.argv[0].split("/")[-1][:-3] + "__" + timestamp

java_performance_improvement_option = "-Xmx16g"
threads_high = subprocess.check_output(["nproc"])
threads_low = str(int(threads_high) / 4)



##################################################
######---------function definitions---------######
##################################################

total_size_of_created_files = 0
total_size_of_removed_intermediates = 0

def print_size(filename):
	size = get_size(filename)
	global total_size_of_created_files
	total_size_of_created_files += size
	name = filename.split("/")[-1]
	prnt("Size of " + name + ": " + str(size/1024/1024) + "G ("+str(size)+")")
def prnt(s):
	print s;
	log.write(s + "\n")
def file_exists(filename): 
	# file name should include directory, or be in working directory
	filename_list = filename.split("/")
	file_path = ""
	for i in range(0,len(filename_list)-1):
		file_path += filename_list[i]+"/"	
	files_in_dir = subprocess.check_output(["ls", file_path]).split()

	filename = filename.split("/")[-1]


	if filename in files_in_dir:
		prnt("File " + filename + " exists. To redo this step, delete the file first.")
		return True
	else:
		prnt("Creating file: " + filename)
		return False
def error_from_tool(tool, return_code):
	command_args = tool.split(" ")

	prnt("Error with " + tool + "! Exited with return code: " + str(return_code))
	if return_code == 1:
		prnt("Please make sure that all index files are present.")
	if "AnalyzeCovariates" in command_args:
		what_next = raw_input("Abort? (enter 'no' to continue) ")
		if what_next == "no":
			prnt("Continuing after Plotting error. Please see tables.")
			return;
	prnt("Aborting.")
	prnt("See logfile: " + log.name)
	sys.exit(1)
def benchmark(first_last="neither"):
	global last_time
	if first_last == "last":
		time_now = datetime.now()
		prnt("Done. Time taken for this step: " + time_diff(time_now,last_time))
		prnt("---")
		prnt("All done. Total time elapsed: " + time_diff(time_now,start_time))
	elif first_last == "first":
		last_time = start_time
		time_now = datetime.now()
        	prnt("Time elapsed since start: " + time_diff(time_now,start_time))
        	last_time = time_now
	elif first_last == "neither":
		time_now = datetime.now()
        	prnt("Done. Time taken for this step: " + time_diff(time_now,last_time))
        	prnt("Time elapsed since start: " + time_diff(time_now,start_time))
        	last_time = time_now
	else:
		prnt("Error in benchmark(): arg = " +first_last)
		sys.exit(1)
def time_diff(a, b):
	if a>b:
		d = a-b
	else:
		d = b-a
	return str(d.days) +" days, "+  str(d.seconds/3600) +":"+ str((d.seconds%3600)/60) + ":" + str((d.seconds%3600)%60) +" (" + str(d.seconds+d.days*86400) + " seconds)"
def rm_intermediates(*file_list):
	global total_size_of_removed_intermediates
	for file in file_list:
		total_size_of_removed_intermediates += get_size(file)
		#subprocess.call(["rm",  file])
		prnt("Not deleting: " + file)
def index_bam(bam_file):
	prnt("Indexing bam file...")
	if not file_exists(bam_file + ".bai"):
		call(["samtools", "index", bam_file])
def call(command, benchmark_arg="neither", stdout="stdout"):
	# prnt(command and save command line version of command )
	#   in case of error message
	prnt("---- Command: ")
	c = ""
	for arg in command:
		c += arg + " "
	if stdout is not "stdout":
		c += " > " + stdout.name
	prnt(c)
	
	# call subprocess, checking for errors 
	if stdout=="stdout":
		return_code = subprocess.call(command, stderr=log) 
	else:
		return_code = subprocess.call(command, stdout=stdout, stderr=log)
	if return_code!= 0:
		error_from_tool(c, return_code)
	
	# benchmark this command
	if benchmark_arg is not "none":
		benchmark(first_last=benchmark_arg)
def get_size(filename):
	return int(subprocess.check_output(["ls", "-s", filename]).split()[0])


##################################################
######-------check pipline arguments--------######
##################################################

###-----
sample_name = reads1.split("/")[-1].split(".")[0]
working_dir = reads1.split(".")[0] + "___" + comment + "/"

# check if "create new folder?" parameter was provided
if len(sys.argv) == 9:
	if sys.argv[8] == "0":
		working_dir = subprocess.check_output(["pwd"]).strip() + "/"
else:
	print "Making new directory: " + working_dir
	subprocess.call(["mkdir", working_dir])
print "\nWorking directory is: " + working_dir

tmp_dir = working_dir + "tmp/"
subprocess.call(["mkdir", tmp_dir])
java_tmp_dir = "-Djava.io.tmpdir=" + tmp_dir
print "Temp directory is: " + tmp_dir


log = open(working_dir + sample_name + "_" + comment + ".log", "w")

prnt("")
prnt("--------------BEGIN-------------")
 
# check and fix pathnames for tools
picard = picard.strip()
snpeff = snpeff.strip()
if picard[-1] != "/":
	picard = picard + "/"
if snpeff[-1] != "/":
	snpeff = snpeff + "/"


# check ref directory
ref_list = ref.split("/")
ref_path = ""
for i in range(0,len(ref_list)-1):
	ref_path += ref_list[i]+"/"

files_in_ref_dir = subprocess.check_output(["ls", ref_path]).split()

# checking if fasta index exists, create if not
if ref_list[-1].split(".")[0] + ".fai" in files_in_ref_dir:
	# fasta index with bad file extension found
	prnt("Fasta index exists: " + ref_list[-1].split(".")[0] + ".fai")
	# if the festa index is ref.fai, bwa won't recognize it. 
	#    it has to be renamed to ref.fa.fai
	prnt("Renaming to: " + ref_list[-1] + ".fa.fai")
    	call(["mv", ref_list[-1] + ".fai", ref_list[-1] + ".fa.fai"], benchmark_arg="none")
elif ref_list[-1] + ".fai" in files_in_ref_dir:
	# correct fasta index found
	prnt("Fasta index exists: " + ref_list[-1] + ".fai")
else:
	# no fasta index found
	prnt("Fasta index does not exist yet. Creating...")
	call(["samtools", "faidx",  ref], benchmark_arg="none")

# checking if fasta dict exists, create if not
if ref_list[-1].split(".")[0] + ".dict" in files_in_ref_dir:
	# fasta dict found
	prnt("Reference sequence dictionary exists: " + ref_list[-1].split(".")[0] + ".dict")
	benchmark(first_last="first")
elif ref_list[-1] + ".dict" in files_in_ref_dir:
	# fasta dict with bad file extension found
	prnt("Reference sequence dictionary exists: " + ref_list[-1] + ".dict")
	# if the festa index is ref.fa.dict, GATK won't recognize it. 
	#    it has to be renamed to ref.dict
	prnt("Renaming to: " + ref_list[-1].split(".")[0] + ".dict")
    	call(["mv", ref_list[-1] + ".dict", ref_list[-1].split(".")[0] + ".dict"],benchmark_arg="first")
else:
	# no fasta dict found
	prnt("Reference sequence dictionary does not exist yet. Creating...")
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", picard+"CreateSequenceDictionary.jar", "REFERENCE="+ref, "OUTPUT="+ref_path+ref_list[-1].split(".")[0] + ".dict"],benchmark_arg="first")


# ref_path is now the path to the directory containing the reference
ref = ""
for i in range(0,len(ref_list[-1].split('.'))-1):
	ref += ref_list[-1].split('.')[i]
ref = ref_path + ref
# ref is now the name of the reference file (without .fa or .fasta)




##################################################
######---------- MAIN: pipeline ------------######
##################################################





# aligning fastq reads with reference
prnt("")
prnt ("-" * 40)
prnt("Aligning reads with reference...")

#reads1_sai = working_dir + reads1.split("/")[-1] + ".sai"
#if not file_exists(reads1_sai):
#	call(["bwa", "aln", "-t"+threads_high,"-f", reads1_sai, ref+".fa", reads1])
#print_size(reads1_sai)
#
#reads2_sai = working_dir + reads2.split("/")[-1] + ".sai"
#if not file_exists(reads2_sai):
#	call(["bwa", "aln", "-t"+threads_high,"-f", reads2_sai, ref+".fa", reads2])
#print_size(reads2_sai)


###-----

# creating sam file
#prnt("")
#prnt("--------------------------------" )
#prnt("Writing alignments to sam file...")
sam_file = working_dir + reads1.split("/")[-1].split(".")[0] + ".sam"
prnt("samfile: " + sam_file)
if not file_exists(sam_file):
	f = open(sam_file, "w")
	#call(["bwa", "sampe", "-f", sam_file, ref+".fa", reads1_sai, reads2_sai, reads1 ,reads2])
	call(["bwa", "mem", "-t"+threads_high, "-M", ref+".fa", reads1, reads2], stdout=f)
	f.close()
#rm_intermediates(reads1_sai, reads2_sai)
print_size(sam_file)

# cleaning up sam file
prnt("")
prnt("--------------------------------" )
prnt("Cleaning sam file...")
previous_sam_file = sam_file
sam_file = previous_sam_file + ".clean.sam"
if not file_exists(sam_file):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", picard + "CleanSam.jar", "INPUT="+previous_sam_file, "OUTPUT="+sam_file])
rm_intermediates(previous_sam_file)
print_size(sam_file)


# converting to  bam file
prnt("")
prnt("--------------------------------" )
prnt("Converting to bam file...")
bam_file = sam_file[:len(sam_file)-3] + "bam"
if not file_exists(bam_file):
	call(["samtools", "view", "-@"+threads_low, "-Sb", "-o", bam_file, sam_file])
rm_intermediates(sam_file)
print_size(bam_file)


# sorting bam file
prnt("")
prnt("--------------------------------" )
prnt("Sorting bam file...")
previous_bam_file = bam_file 
bam_file = previous_bam_file + ".sort"
if not file_exists(bam_file + ".bam"):
	call(["samtools", "sort", "-@"+threads_low, "-m4G", previous_bam_file, bam_file])
rm_intermediates(previous_bam_file)
bam_file = bam_file + ".bam"
print_size(bam_file)

# indexing bam file
index_bam(bam_file)


###-----

# adding readgroups
prnt("")
prnt("--------------------------------" )
prnt("Adding read groups...")
previous_bam_file = bam_file
bam_file = previous_bam_file + ".rg.bam"
if not file_exists(bam_file):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", picard+"AddOrReplaceReadGroups.jar", "I="+previous_bam_file, "O="+bam_file, "RGLB=lib"+sample_name, "RGPL=illumina", "RGPU=barcode", "RGSM="+sample_name])
rm_intermediates(previous_bam_file)
print_size(bam_file)

# indexing bam file
index_bam(bam_file)

# removing duplicates
prnt("")
prnt("--------------------------------" )
prnt("Removing duplicates...")
previous_bam_file = bam_file
bam_file = previous_bam_file + ".rmdup.bam"
if not file_exists(bam_file):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", picard+"MarkDuplicates.jar", "I="+previous_bam_file, "O="+bam_file, "METRICS_FILE="+bam_file+".metrics_file", "REMOVE_DUPLICATES=true", "ASSUME_SORTED=true", "VALIDATION_STRINGENCY=LENIENT"])
rm_intermediates(previous_bam_file)
print_size(bam_file)

# indexing bam file
index_bam(bam_file)


###-----

# indel realigner target creator
prnt("")
prnt("--------------------------------" )
prnt("Finding indels...")
realignment_intervals = bam_file + ".intervals"
if not file_exists(realignment_intervals):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", gatk, "-I", bam_file, "-R", ref+".fa", "-T", "RealignerTargetCreator", "-o", realignment_intervals])
print_size(realignment_intervals)

# indel realigner
prnt("")
prnt("--------------------------------" )
prnt("Realigning around indels...")
previous_bam_file = bam_file
bam_file = previous_bam_file + ".realigned.bam"
if not file_exists(bam_file):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", gatk, "-I",previous_bam_file, "-R", ref+".fa", "-T", "IndelRealigner", "-o", bam_file, "-targetIntervals", realignment_intervals])
rm_intermediates(previous_bam_file)
print_size(bam_file)


# indexing bam file
index_bam(bam_file)

################################################
####--------END THIS PIPELINE VERSION--------###
#prnt("This is where the pipeline ends for now. Duplicate Marking and Indel Realignment are done; Base Score Recalibration wasn't done yet.")
#sys.exit(0)
###############################################


###-----

# base recal step 1
prnt("")
prnt("--------------------------------" )
prnt("Recalibrating base quality scores (Step 1)...")
recal_table1 = bam_file + ".firstpass.table"
if not file_exists(recal_table1):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", gatk, "-I", bam_file, "-R", ref+".fa", "-T", "BaseRecalibrator", "-o", recal_table1, "-knownSites", knownsites])
print_size(recal_table1)

prnt("")
prnt("--------------------------------" )
prnt("Recalibrating base quality scores (Step 2)...")
recal_table2 = bam_file + ".secondpass.table"
if not file_exists(recal_table2):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", gatk, "-I", bam_file, "-R", ref+".fa", "-T", "BaseRecalibrator", "-o", recal_table2, "-knownSites", knownsites, "-BQSR", recal_table1])
print_size(recal_table2)

prnt("")
prnt("--------------------------------" )
prnt("Plotting recalibration data...")
recal_plots = bam_file + ".recal.pdf"
if not file_exists(recal_plots):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", gatk, "-R", ref+".fa", "-T", "AnalyzeCovariates", "-before", recal_table1, "-after", recal_table2, "-plots", recal_plots])
#print_size(recal_plots)

#base recal step 1
prnt("")
prnt("--------------------------------" )
prnt("Writing recalibrated alignments...")
previous_bam_file = bam_file
bam_file = previous_bam_file + ".recal.bam"
if not file_exists(bam_file):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", gatk, "-I", previous_bam_file, "-R", ref+".fa", "-T", "PrintReads", "-o", bam_file, "--BQSR", recal_table2])
print_size(bam_file)

###-----

# calling variants
prnt("")
prnt("--------------------------------" )
prnt("Calling SNPs...")
raw_vcf = bam_file + ".raw.vcf"
if not file_exists(raw_vcf):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", gatk, "-I", bam_file, "-R", ref+".fa", "-T", "UnifiedGenotyper", "-o", raw_vcf, "--dbsnp", knownsites])
print_size(raw_vcf)


###-----

# filtering snps
prnt("")
prnt("--------------------------------" )
prnt("Filtering SNPs...")
filtered_vcf = bam_file + ".filtered.vcf"
if not file_exists(filtered_vcf):
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar",  gatk,  "--variant",  raw_vcf , "-R",  ref+".fa" , "-T", "VariantFiltration", "-o" , filtered_vcf])
rm_intermediates(raw_vcf)
print_size(filtered_vcf)


###-----

# annotation
prnt("")
prnt("--------------------------------" )
prnt("Annotating SNPs...")
annotated_vcf = bam_file + ".filtered.annotated.vcf"
if not file_exists(annotated_vcf):
	f = open(annotated_vcf, "w")
	call(["java", java_performance_improvement_option, java_tmp_dir, "-jar", snpeff+"snpEff.jar", "-v", "-c", snpeff+"snpEff.config", "GRCh37.70", "-s",  bam_file+".html", filtered_vcf], benchmark_arg="last", stdout=f)
	f.close()
rm_intermediates(filtered_vcf)
print_size(annotated_vcf)

prnt("")
prnt("Total size of files created: " + str(total_size_of_created_files))
prnt("Total size of intermediates removed: " + str(total_size_of_removed_intermediates))

prnt("--------------------------------" )
prnt("Log file: " + str(log.name))
prnt("")
prnt("--------------END---------------")
prnt("")
log.close()




