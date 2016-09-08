#!/usr/bin/env python
#	ToDo List:
#			- Fix gff opening
#			- bam processing


from Chromosom import Chromosom as Chrom
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import pandas as pd
import sys
import numpy as np
import threading
import timeit
import os

'''
chromdict: {Chr1 : class(Chromosom), }
int(Chromlength)
vcf = [SNPpos, SNPpos,...]
gff = [[Start, Stop, Info],...]
nlist = [[Start, Stop],...]
desertlist = [Start, Stop, Length, [SNPpos,...],...]
genelist = [Start, Stop, Info, Sequence, [SNPpos,...],...]
'''


def print_help():
	helpline = """
USAGE:

	Oasis -vcf <vcf-path> -fasta <fasta-path> -gff <gff-path> -depth <depth-path> -contig <agp-file> <options>

MANDATORY OPTIONS:

	-f:	Number of SNPs allowed in 1 kB (default: 0)
	-s:	Size of the initial desert seed (in bp) without SNPs (default: 1000)
	-d:	Minimal size of a SNP desert to be considered (default: seed size (1000))

ADDITIONAL OPTIONS:

	-depth:	Path to the file containing the sequence coverage. If it is provided deserts are filtered with the settings of -c and -m
	-agp:	Path to the file containing the positions of the individual contigs. The deserts will be cut to fit these contigs
	-out:	Output path (default: location of the vcf file)
	-c:	Minimal coverage for a SNP desert (default: 15)
	-m:	Maximum percentage of Bases without minimal coverage in a SNP desert (default: 0.0)
	-r:	Maximum percentage of 'n'-characters in a seed (default: 0.0)
	-n:	Maximum percentage of 'n'-characters in a SNP desert (default: 0.0)
	-g:	Removes all genes with a SNP
	-p:	Plots additional graphics
	-h:	Prints this helpline
"""
	print helpline
	return None

def open_file(filepath):
	text = open(filepath, 'r')
	lines = text.readlines()
	text.close()
	filelist = list(lines)
	return filelist

def write_file(filepath, filelist):
	output = open(filepath, 'w')
	output.writelines(filelist)
        output.close
	return None

def open_fasta(fastapath, chromdict):
	fastalist = ''.join(open_file(fastapath)).split('>')
	for Chromseq in fastalist:
		if len(Chromseq) > 0:
			Chromosom = Chromseq.split('\n')[0].split('\t')[0].split(' ')[0]
			Sequence = ''.join(Chromseq.split('\n')[1:])
			chromdict[Chromosom] = Chrom(Chromosom)
			chromdict[Chromosom].set_sequence(Sequence)
	return chromdict

def open_vcf(vcfpath, chromdict, helpdict):
	vcflist = open_file(vcfpath)
	for line in vcflist:
		if line[0] != '#':
			helpdict[line.split('\t')[0]].append(int(line.split('\t')[1]))
	for Chromosom in helpdict:
		chromdict[Chromosom].set_SNPlist(helpdict[Chromosom])
		helpdict[Chromosom] = []
	return None

def open_gff(gffpath, chromdict, helpdict):
	gfflist = open_file(gffpath)
	for locus in gfflist[:-1]:
		if 'mRNA\t' in locus:
			splitline = locus.split('\t')
			Chromosom = splitline[0]
			strand = splitline[6]
			information = splitline[8]
			start = int(splitline[3])-1
			stop = int(splitline[4])-1
			helpdict[Chromosom].append([start, stop, strand, information])
	for Chromosom in helpdict:
		chromdict[Chromosom].set_genelist(helpdict[Chromosom])
		helpdict[Chromosom] = []
	return None

def open_depth(depthpath, mincoverage, chromdict, helpdict):
	for Chromosom in helpdict:
		helpdict[Chromosom] = [[],0]
	infile = open(depthpath, 'r' )
	reader = csv.reader(infile, delimiter = '\t')
	for line in reader:
		if len(line[0].split(':')) > 1:
			helpdict[line[0].split(':')[0]][1] += int(line[1])
			if int(line[1]) < mincoverage:
				helpdict[line[0].split(':')[0]][0].append(int(line[0].split(':')[1])-1)
	for Chromosom in helpdict:
		chromdict[Chromosom].set_coveragelist(helpdict[Chromosom][0])
		chromdict[Chromosom].set_coverage(helpdict[Chromosom][1])
		helpdict[Chromosom] = []
	return None

def open_bam(bampath, mincoverage, chromdict, helpdict):
	return None

def open_agp(agppath, chromdict, helpdict):
	for Chromosom in helpdict:
		helpdict[Chromosom] = []
	agplist = open_file(agppath)
	for line in agplist:
		splitline = line.split('\t')
		if splitline[3] != 'N':
			start = int(splitline[1])
			stop = int(splitline[2])
			helpdict[splitline[0]].append([start, stop])
	for Chromosom in helpdict:
		chromdict[Chromosom].set_contiglist(helpdict[Chromosom])
	return None
	

def multithread_chromosom(chromdict, sensitivity, seedsize, minsize, desertremove, desertcutoff, seedremove, seedcutoff, covcutoff, desertfitting):
	threadlist = []
	for Chromosom in chromdict:
	 	t = threading.Thread(target=chromdict[Chromosom].process_data, args=(sensitivity, seedsize, minsize, desertremove, desertcutoff, seedremove, seedcutoff, covcutoff, desertfitting,))
  		threadlist.append(t)
	for thread in threadlist:
		thread.start()
	for thread in threadlist:
		thread.join()
	return None

def write_output_files(chromdict, outpath):
	desertlines = ['##gff-version 3\n']; genelines = ['##gff-version 3\n']; sizelist = []
	desertpath = outpath+'_deserts.gff3'
	genepath = outpath+'_genes.gff3'
	histogrampath = outpath+'_hist.png'
	chromgraphpath = outpath+'_chrom.png'
	numberofdeserts = []; desertpercentage = []; chromnames = []
	for Chromosom in chromdict:
		lines, sizes, desertcov =  make_desert_lines(chromdict[Chromosom])
		desertlines += lines
		sizelist += sizes
		lines =  make_gene_lines(chromdict[Chromosom])
		genelines += lines
		chromnames.append(chromdict[Chromosom].get_name())
		numberofdeserts.append(len(chromdict[Chromosom].get_desertlist()))
		desertpercentage.append(desertcov)
		if plots == True:
			chromdict[Chromosom].make_desertplot(outpath)
	make_histogram(sizelist, histogrampath)
	make_chromosomgraph(numberofdeserts, desertpercentage, chromnames, chromgraphpath)
	write_file(desertpath, desertlines)
	write_file(genepath, genelines)
	return None

def make_histogram(sizelist, histogrampath):
	sizehistogram, bins = np.histogram(sizelist, bins=80)
	width = 0.7 * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	plt.title('sizes of SNP deserts')
	plt.ylabel('# of deserts')
	plt.xlabel('sizes [kB]')
	plt.bar(center, sizehistogram, align='center', width=width, color='b', log=True)
	plt.savefig(histogrampath, dpi=150)
	return None

def make_chromosomgraph(numberofdeserts, desertpercentage, chromnames, chromgraphpath):
	fig = plt.figure(); ax = fig.add_subplot(111)
	N = len(chromdict); ind = np.arange(N); width = 0.35 

	rects1 = ax.bar(ind, numberofdeserts, width, color='black')
	rects2 = ax.bar(ind+width, desertpercentage, width, color='red')

	# axes and labels
	ax.set_xlim(-width,len(ind)+width)
	#ax.set_ylim(0,45)
	ax.set_title('# of SNP deserts per chromosom & percentage of covered sequence')
	xTickMarks = chromnames
	ax.set_xticks(ind+width)
	xtickNames = ax.set_xticklabels(xTickMarks)
	plt.setp(xtickNames, rotation=45, fontsize=10)

	## add a legend
	ax.legend( (rects1[0], rects2[0]), ('# of SNP deserts', '% of covered Sequence') )

	plt.savefig(chromgraphpath, dpi=150)
	return None

def make_desert_lines(Chromosom):
	linelist = ['##sequence-region '+Chromosom.get_name()+' 1 '+str(Chromosom.get_size())+'\n']; desertlength = 0; sizelist = []; i = 1
	for desert in Chromosom.get_desertlist():
		line = '\t'.join([Chromosom.get_name(),'Oasis','Region',str(desert[0]+1),str(desert[1]+1),'.','.','.','Note='+str(len(desert[3]))+'\n'])
		desertlength += desert[2]
		sizelist.append(desert[2]/1000.0)
		linelist.append(line)
	desertpercentage = desertlength/float(Chromosom.get_size())*100
	return linelist, sizelist, desertpercentage
	
def make_gene_lines(Chromosom):
	linelist = ['##sequence-region '+Chromosom.get_name()+' 1 '+str(Chromosom.get_size())+'\n']; i = 1
	for gene in Chromosom.get_desertgenes():
		line = '\t'.join([Chromosom.get_name(),'Oasis','gene',str(gene[0]+1),str(gene[1]+1),'.',gene[2],'.',gene[3]])
		if generemove == True:
			if len(gene) == 4:
				linelist.append(line)
		else:
			linelist.append(line)
		i+=1
	return linelist
		
def get_SNP_spacing(chromdict, helpdict, outpath):
	histogrampath = outpath+'_space.png'; spacings = []
	for chrom in chromdict:
		spacings.append(chromdict[chrom].get_SNP_spacing())
	fig, ax1 = plt.subplots(figsize=(10,6))
	fig.canvas.set_window_title('A Boxplot Example')
	plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

	bp = plt.boxplot(spacings, whis=1000000)
	plt.yscale('log')
	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='red', marker='o')
	plt.savefig(histogrampath, dpi=150)
	spacelist = []; SNPlist = []; listpath = outpath+'_spacelist.lst'
#	for part in spacings:
#		for space in part:
#			spacelist.append(str(space))
#	write_file(listpath, '\n'.join(spacelist))
	return None
	
'''
	spacehistogram, bins = np.histogram(spacelist, bins=100)
	width = 0.7 * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	plt.title('base pairs between SNPs')
	plt.ylabel('# of SNPs')
	plt.xlabel('sizes [kB]')
	plt.bar(center, spacehistogram, align='center', width=width, color='b', log=True)
	plt.savefig(histogrampath, dpi=150)
	strlist = [str(space)+'\n' for space in spacelist]
	write_file(listpath, strlist)
	return None
	'''	


arguments = sys.argv[1:]
covcutoff = 0
desertcutoff = 0
seedcutoff = 0
sensitivity = 0
seedsize = 1000
minsize = 1000
mincoverage = 15
desertremove = False
seedremove = False
generemove = False
desertfitting = False
plots = False

if len(sys.argv) <= 1:
	print_help()
else:
	if '-vcf' in sys.argv:
		if os.path.isfile(sys.argv[sys.argv.index('-vcf')+1]):
			vcfpath = sys.argv[sys.argv.index('-vcf')+1]
		else:
			sys.exit("ERROR: no vcf-file provided!\n")
	else:
		sys.exit("ERROR: no vcf-file provided!\n")
	if '-gff' in sys.argv:
		if os.path.isfile(sys.argv[sys.argv.index('-gff')+1]):
			gffpath = sys.argv[sys.argv.index('-gff')+1]
		else:
			sys.exit("ERROR: no gff-file provided!\n")
	else:
		sys.exit("ERROR: no gff-file provided!\n")
	if '-fasta' in sys.argv:
		if os.path.isfile(sys.argv[sys.argv.index('-fasta')+1]):
			fastapath = sys.argv[sys.argv.index('-fasta')+1]
		else:
			sys.exit("ERROR: no fasta-file provided!\n")
	else:
		sys.exit("ERROR: no fasta-file provided!\n")
	if '-depth' in sys.argv:
		if os.path.isfile(sys.argv[sys.argv.index('-depth')+1]):
			covpath = sys.argv[sys.argv.index('-depth')+1]
	else:
		print "WARNING: no depth-file provided!\nNo coverage will be calculated"
		covpath = False
	if '-contig' in sys.argv:
		if os.path.isfile(sys.argv[sys.argv.index('-contig')+1]):
			agppath = sys.argv[sys.argv.index('-contig')+1]
			desertfitting = True
	else:
		print "WARNING: no agp-file provided!\nNo contigs will be separated"
		agppath = False
	if '-f' in sys.argv:
		try:
			sensitivity = int(sys.argv[sys.argv.index('-f')+1])
		except:
			sys.exit("ERROR: wrong value given for seed-SNP-frequency!\n")
	if '-s' in sys.argv:
		try:
			seedsize = int(sys.argv[sys.argv.index('-s')+1])
		except:
			sys.exit("ERROR: wrong value given for seedsize!\n")
	if '-c' in sys.argv:
		try:
			mincoverage = int(sys.argv[sys.argv.index('-c')+1])
		except:
			sys.exit("ERROR: wrong value given for minimal coverage!\n")
	if '-d' in sys.argv:
		try:
			minsize = int(sys.argv[sys.argv.index('-d')+1])
		except:
			minsize = seedsize
	if '-m' in sys.argv:
		try:
			covcutoff = float(sys.argv[sys.argv.index('-m')+1])
		except:
			sys.exit("ERROR: wrong value given for coverage cutoff!\n")
	if '-n' in sys.argv:
		try:
			desertcutoff = float(sys.argv[sys.argv.index('-n')+1])
			desertremove = True
		except:
			sys.exit("ERROR: wrong value given for desert 'n' character cutoff!\n")
	if '-r' in sys.argv:
		try:
			seedcutoff = float(sys.argv[sys.argv.index('-r')+1])
			seedremove = True
		except:
			sys.exit("ERROR: wrong value given for seed 'n' character cutoff!\n")
	if not '-d' in sys.argv:
		minsize = seedsize
	if '-h' in sys.argv:
		print_help()
	if '-g' in sys.argv:
		generemove = True
	if '-p' in sys.argv:
		plots = True
	if '-out' in sys.argv:
		try:
			outpath = sys.argv[sys.argv.index('-out')+1]
		except:
			sys.exit("ERROR: no valid output path provided!\n")
	if not '-out' in sys.argv:
		try:
			outpath = '.'.join(vcfpath.split('.')[0:-1])+'_oasis_f'+str(sensitivity)
			print 'output directed to: '+outpath
		except:
			sys.exit("ERROR: no valid output path provided!\n")

	chromdict = {}; helpdict = {}
	print 'processing fasta file'
	chromdict = open_fasta(fastapath, chromdict)
	for Chromosom in chromdict:
		helpdict[Chromosom] = []
	print 'processing vcf file'
	open_vcf(vcfpath, chromdict, helpdict)
	print 'processing gff file'
	open_gff(gffpath, chromdict, helpdict)
	if covpath != False:
		print 'processing depth file'
		if '.bam' in covpath.lower():
			open_bam(covpath, mincoverage, chromdict, helpdict)
		else:
			open_depth(covpath, mincoverage, chromdict, helpdict)
	if agppath != False:
		print 'processing agp-file'
		open_agp(agppath, chromdict, helpdict)
	print 'multithreading chromosoms'
	multithread_chromosom(chromdict, sensitivity, seedsize, minsize, desertremove, desertcutoff, seedremove, seedcutoff, covcutoff, desertfitting)
	print 'writing files'
	write_output_files(chromdict, outpath)
	
	#get_SNP_spacing(chromdict, helpdict, outpath)

#start = timeit.default_timer()
#stop = timeit.default_timer()
#print stop - start 

