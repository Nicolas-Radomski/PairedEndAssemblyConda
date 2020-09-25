#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

#### author: Nicolas Radomski ####
# version for conda in a cluster: conda enviroment PairedEndAssembly
# run BBnorn (step 1_normalization), Trimmomatic (step 2_trimming), Spades (step 3_assembly), prokka (4_annotation) and Quast (5_quality) 
# based on paired-end reads from a single genomic sample
# the present main script PairedEndAssemblyConda.py and module genomic.py (version 20200923) were prepared and tested with Python and Conda packages below (Septembre 2020)
# Name                    Version                   Build  Channel
# bbmap                     38.84                h516909a_0    bioconda
# biopython                 1.78             py37h8f50634_0    conda-forge
# python                    3.7.8           h6f2ec95_1_cpython    conda-forge
# quast                     5.0.2           py37pl526hb5aa323_2    bioconda
# spades                    3.13.1                        0    bioconda
# trimmomatic               0.39                          1    bioconda
# prokka                    1.14.6                  pl526_0    bioconda
# the module genomic.py has to be with the present main script PairedEndAssemblyConda.py to lunch it
# Execute /global/conda/envs/PairedEndAssembly/bin/bbnorm.sh in=/global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/ERR3997398_R1.fastq.gz in2=/global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/ERR3997398_R2.fastq.gz out=/global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/1_normalization/ERR3997398_R1_N.fastq.gz out2=/global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/1_normalization/ERR3997398_R2_N.fastq.gz target=100 threads=48
# Execute /global/conda/envs/PairedEndAssembly/bin/trimmomatic PE -threads 48 -phred33 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/1_normalization/ERR3997398_R1_N.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/1_normalization/ERR3997398_R2_N.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/2_trimming/ERR3997398_R1_P.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/2_trimming/ERR3997398_R1_UP.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/2_trimming/ERR3997398_R2_P.fastq.gz /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/2_trimming/ERR3997398_R2_UP.fastq.gz ILLUMINACLIP:/global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# Execute /global/conda/envs/PairedEndAssembly/bin/spades.py -1 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/2_trimming/ERR3997398_R1_P.fastq.gz -2 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/2_trimming/ERR3997398_R2_P.fastq.gz -t 48 --phred-offset 33 -o /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/3_assembly/
# Execute /global/conda/envs/PairedEndAssembly/bin/prokka --outdir /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/4_annotation/ --force --prefix ERR3997398 --addgenes --cpus 48 --proteins /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109_prot.fasta /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/3_assembly/ERR3997398_contigs.header.single.number.parsed.fasta
# Execute /global/conda/envs/PairedEndAssembly/bin/quast.py -o /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/5_quality/ -r /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.fasta -g /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.gff3 -m 500 -t 48 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/ERR3997398/3_assembly/ERR3997398_contigs.header.single.number.parsed.fasta

'''
#### exemple of Bash command (sbatch_PairedEndAssemblyConda.sh) ####
#!/bin/bash
#SBATCH -p Research
#SBATCH -o %x.%N.%j.out
#SBATCH -e %x.%N.%j.err
#SBATCH --cpus-per-task=48
#SBATCH --job-name=test-20200922
source /global/conda/bin/activate;conda activate PairedEndAssembly; \
python /global/bio/projets/GAMeR/Nicolas-Radomski/Python/PairedEndAssemblyConda.py \
 -t 48 \
 -n 100 \
 -l 500 \
 -R1 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/ERR3997398_R1.fastq.gz \
 -R2 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/ERR3997398_R2.fastq.gz \
 -adap /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/NexteraPE-PE.fa \
 -prot /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109_prot.fasta \
 -nucl /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.fasta \
 -coor /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.gff3 \
 -norm /global/conda/envs/PairedEndAssembly/bin/bbnorm.sh \
 -trim /global/conda/envs/PairedEndAssembly/bin/trimmomatic \
 -denovo /global/conda/envs/PairedEndAssembly/bin/spades.py \
 -annot /global/conda/envs/PairedEndAssembly/bin/prokka \
 -qual /global/conda/envs/PairedEndAssembly/bin/quast.py
#### exemple of Bash command execution ####
sbatch -p Research --cpus-per-task=48 -e %x.%N.%j.err -o %x.%N.%j.out sbatch_PairedEndAssemblyConda.sh
'''

import os, sys
import argparse
import genomic

# parse arguments
def get_parser():
	
	# function asking arguments
	parser = argparse.ArgumentParser(description="perform read normalization, read trimming, de novo assembly, contig annotation and contig quality from forward (ID_R1.fastq.gz) and reverse (ID_R2.fastq.gz) paired-end reads (offset of reads must be Phred 33 and only uppercase or lowercase alphanumeric characters are expected in the ID sample)")

	# setting of arguments
	parser.add_argument('-t', action="store", dest='threads',
					type=int, required=True, 
					help='number of threads (REQUIRED)')

	parser.add_argument('-n', action="store", dest='depth',
					type=int, required=True, 
					help='read normalization depth (REQUIRED)')
					
	parser.add_argument('-l', action="store", dest='length',
					type=int, required=True, 
					help='minimal length of contigs (REQUIRED)')

	parser.add_argument('-R1', action="store", dest='forward',
					type=str, required=True, 
					help='path to the forward read (REQUIRED)')

	parser.add_argument('-R2', action="store", dest='reverse',
					type=str, required=True, 
					help='path to the reverse read (REQUIRED)')

	parser.add_argument('-adap', action="store", dest='adapters',
					type=str, required=True, 
					help='path to adapters for trimming (REQUIRED)')

	parser.add_argument('-prot', action="store", dest='proteins',
					type=str, required=True, 
					help='path to reference protein fasta for annotation (REQUIRED)')

	parser.add_argument('-nucl', action="store", dest='nucleotides',
					type=str, required=True, 
					help='path to reference nucleotide fasta for quality (REQUIRED)')

	parser.add_argument('-coor', action="store", dest='coordinates',
					type=str, required=True, 
					help='path to reference coordinate gff for quality (REQUIRED)')

	parser.add_argument('-norm', action="store", dest='bbnorm',
					type=str, required=True, 
					help='path to BBnorm (REQUIRED)')					

	parser.add_argument('-trim', action="store", dest='trimmomatic',
					type=str, required=True, 
					help='path to Trimmomatic (REQUIRED)')

	parser.add_argument('-denovo', action="store", dest='spades',
					type=str, required=True, 
					help='path of SPAdes (REQUIRED)')

	parser.add_argument('-annot', action="store", dest='prokka',
					type=str, required=True, 
					help='path to Prokka (REQUIRED)')

	parser.add_argument('-qual', action="store", dest='quast',
					type=str, required=True, 
					help='path to Quast (REQUIRED)')

	return parser

# ask threads, depth of read normalization, contig length, R1 path, R2 path, adapter path, 
# path of reference proteins, path of reference nucleotids, path of reference coordinates,
# BBnorm path, Trimmomatic path, SPAdes path, Prokka path, Quast path,
# then return directories 1_normalization, 2_trimming, 3_assembly, 4_annotation and 5_quality
def main():
	
	# get parser object
	parser=get_parser()

	# print parser.help if there are no arguments in the command
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)

	# extract arguments from parser
	Arguments=parser.parse_args()

	# define variables of the arguments	
	t=Arguments.threads
	n=Arguments.depth
	l=Arguments.length
	R1=Arguments.forward
	R2=Arguments.reverse
	AD=Arguments.adapters
	RP=Arguments.proteins
	RN=Arguments.nucleotides
	RC=Arguments.coordinates
	BB=Arguments.bbnorm
	TR=Arguments.trimmomatic
	SP=Arguments.spades
	PR=Arguments.prokka
	QU=Arguments.quast

	# get and print working directory (wd)
	wd = os.getcwd()
	print ("#### The current working directory is %s" % wd)

	# extract sample name and ID
	R1name = os.path.basename(R1)
	print("#### The forward read is named %s" %R1name)
	R2name = os.path.basename(R2)
	print("#### The reverse read is named %s" %R2name)
	R1id = R1name.replace('_R1.fastq.gz','')
	print("#### The forward ID is %s" %R1id)
	R2id = R2name.replace('_R2.fastq.gz','')
	print("#### The reverse ID is %s" %R2id)
	
	# check IDs identity, then prepare output directories and run BBnorm
	if R1id != R2id:
		print("#### The forward (ID_R1.fastq.gz) and reverse (ID_R2.fastq.gz) reads do not match")
	else:
		print("#### The forward (ID_R1.fastq.gz) and reverse (ID_R2.fastq.gz) reads match")
		# create an output directory called ID sample if it does not exist with the function nodir_makedir_warning of the module genomic.py
		genomic.nodir_makedir_warning(directory = R1id)
		# create output directories if they do not exist with the function nodir_makedir_warning of the module genomic.py
		genomic.nodir_makedir_warning(directory = R1id + '/' + '1_normalization')
		genomic.nodir_makedir_warning(directory = R1id + '/' + '2_trimming')
		genomic.nodir_makedir_warning(directory = R1id + '/' + '3_assembly')
		genomic.nodir_makedir_warning(directory = R1id + '/' + '4_annotation')
		genomic.nodir_makedir_warning(directory = R1id + '/' + '5_quality')
		# prepare path of output directories
		normalizationoutput = wd + '/' + R1id + '/' + '1_normalization' + '/'
		trimmingoutput = wd + '/' + R1id + '/' + '2_trimming' + '/'
		assemblyoutput = wd + '/' + R1id + '/' + '3_assembly' + '/'
		annotationoutput = wd + '/' + R1id + '/' + '4_annotation' + '/'
		qualityoutput = wd + '/' + R1id + '/' + '5_quality' + '/'
		# prepare and run BBnorm
		R1N = normalizationoutput + R1id + '_R1_N.fastq.gz'
		R2N = normalizationoutput + R2id + '_R2_N.fastq.gz'
		cmdBBnorm = BB + ' in=' + R1 + ' in2=' + R2 + ' out=' + R1N + ' out2=' + R2N + ' target=' + str(n) + ' threads=' + str(t)
		print("#### Execute %s" %cmdBBnorm)
		os.system(cmdBBnorm)

	# check the step 1_normalization, then return warning or successfull message with the function emptydir_warning_success of the module genomics.py
	genomic.emptydir_warning_success(directory = R1id + '/' + '1_normalization')

	# prepare and run Trimmomatic
	R1P = trimmingoutput + R1id + '_R1_P.fastq.gz'
	R1UP  = trimmingoutput + R1id + '_R1_UP.fastq.gz'
	R2P = trimmingoutput + R2id + '_R2_P.fastq.gz'
	R2UP = trimmingoutput + R2id + '_R2_UP.fastq.gz'
	cmdTrimmomatic = TR + ' PE ' + '-threads ' + str(t) + ' -phred33 ' + R1N + ' ' + R2N + ' ' + R1P + ' ' + R1UP + ' ' + R2P + ' ' + R2UP + ' ILLUMINACLIP:' + AD + ':2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
	print("#### Execute %s" %cmdTrimmomatic)
	os.system(cmdTrimmomatic)

	# check the step 2_trimming, then return warning or successfull message with the function emptydir_warning_success of the module genomics.py
	genomic.emptydir_warning_success(directory = R1id + '/' + '2_trimming')
	
	# prepare and run SPAdes
	R1P = trimmingoutput + R1id + '_R1_P.fastq.gz'
	R2P = trimmingoutput + R2id + '_R2_P.fastq.gz'
	SPAdescmd = SP + ' -1 ' + R1P + ' -2 ' + R2P + ' -t ' + str(t) + ' --phred-offset 33' + ' -o ' + assemblyoutput
	print("#### Execute %s" %SPAdescmd)
	os.system(SPAdescmd)
		
	# check the step 3_assembly, then return warning or successfull message with the function emptydir_warning_success of the module genomics.py
	genomic.emptydir_warning_success(directory = R1id + '/' + '3_assembly')

	# add prefix to files produced during the step 3_assembly with the function prefix_files_dir of the module genomics.py
	genomic.prefix_files_dir(directory = R1id + '/' + '3_assembly', prefix = R1id)

	# replace header of the contigs.fasta file by ID the sample with the function string_in_fasta of the module genomics.py
	assemblycontigs = assemblyoutput + R1id + '_' + 'contigs.fasta'
	headercontigs = assemblyoutput + R1id + '_' + 'contigs.header.fasta'
	genomic.string_in_fasta(infilepath = assemblycontigs, outfilepath = headercontigs, ID = R1id)

	# pass the contigs.fasta file from multi to single lines with the function multi_single_fasta of the module genomics.py
	singlelinecontigs = assemblyoutput + R1id + '_' + 'contigs.header.single.fasta'
	genomic.multi_single_fasta(infilepath = headercontigs, outfilepath = singlelinecontigs)

	# add incrementing number at the end of headers of the contigs.fasta file with the function number_in_fasta of the module genomics.py
	numbercontigs = assemblyoutput + R1id + '_' + 'contigs.header.single.number.fasta'
	genomic.number_in_fasta(infilepath = singlelinecontigs, outfilepath = numbercontigs)

	# pass the contigs.fasta file from single to multi lines with the function single_multi_fasta of the module genomics.py
	multilinecontigs = assemblyoutput + R1id + '_' + 'contigs.header.multi.fasta'
	genomic.single_multi_fasta(infilepath = numbercontigs, outfilepath = multilinecontigs)

	# remove small contigs from the single line fasta file with the function parse_contigs_fasta of the module genomics.py
	parsedsinglelinecontigs = assemblyoutput + R1id + '_' + 'contigs.header.single.number.parsed.fasta'
	genomic.parse_contigs_fasta(infilepath = numbercontigs, outfilepath = parsedsinglelinecontigs, length = int(l))

	# pass the contigs.fasta file from single to multi lines with the function single_multi_fasta of the module genomics.py
	parsedmultilinecontigs = assemblyoutput + R1id + '_' + 'contigs.header.multi.number.parsed.fasta'
	genomic.single_multi_fasta(infilepath = parsedsinglelinecontigs, outfilepath = parsedmultilinecontigs)

	# prepare and run Prokka
	prokkainput =  parsedsinglelinecontigs
	Prokkacmd = PR + ' --outdir ' + annotationoutput + ' --force' + ' --prefix ' + R1id + ' --addgenes' + ' --cpus ' + str(t) + ' --proteins ' + RP + ' ' + prokkainput
	print("#### Execute %s" %Prokkacmd)
	os.system(Prokkacmd)

	# check the step 4_annotation, then return warning or successfull message with the function emptydir_warning_success of the module genomic.py
	genomic.emptydir_warning_success(directory = R1id + '/' + '4_annotation')

	# prepare and run Quast
	quastinput =  parsedsinglelinecontigs
	Quastcmd = QU + ' -o ' + qualityoutput + ' -r ' + RN + ' -g ' + RC + ' -m ' + str(l) + ' -t ' + str(t) + ' ' + quastinput
	print("#### Execute %s" %Quastcmd)
	os.system(Quastcmd)

	# add prefix to files produced during the step 5_quality with the function prefix_files_dir of the module genomics.py
	genomic.prefix_files_dir(directory = R1id + '/' + '5_quality', prefix = R1id)

	# check the step 5_quality, then return warning or successfull message with the function emptydir_warning_success of the module genomic.py
	genomic.emptydir_warning_success(directory = R1id + '/' + '5_quality')

	# congtratulate users with the function congratulation of the module genomic.py
	genomic.congratulation()

# driver code: if the code above is a scrypt, call  main() function, rather than to considere it as a module
if __name__ == "__main__":
	# calling main() function
	main()
