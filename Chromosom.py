#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class Chromosom:

	def __init__(self, name):
		self.name = name
		self.size = 0
		self.nfrequency = 0.0
		self.SNPfrequency = 0.0
		self.coverage = 0.0
		self.sequence = ''
		self.SNPlist = [] 	#SNP Position
		self.nlist = []		#[start, stop]
		self.genelist = []	#[start, stop, strand, information]
		self.desertlist = []	#[Start, Stop, Length, [SNPpos,...]]
		self.desertgenes = []	#[Start, Stop, strand, information, Sequence, [SNPpos,...]]
		self.coveragelist = []	#
		self.contiglist = []
		self.desertplot = None

	def get_name(self):
		return self.name

	def set_size(self, size):
		self.size = size

	def get_size(self):
		return self.size

	def set_nfrequency(self, frequency):
		self.nfrequency = frequency

	def get_nfrequency(self):
		return self.nfrequency

	def set_SNPfrequency(self, frequency):
		self.SNPfrequency = frequency

	def get_SNPfrequency(self):
		return self.SNPfrequency

	def set_coverage(self, coverage):
		self.coverage = float(coverage)/self.size

	def get_coverage(self):
		return self.coverage

	def set_sequence(self, sequence):
		self.set_size(len(sequence))
		self.make_nlist(sequence)
		

	def get_sequence(self):
		return self.sequence

	def set_SNPlist(self, SNPlist):
		self.SNPlist = sorted(SNPlist)
		if self.size != 0:
			self.SNPfrequency = len(self.SNPlist)/float(self.size)

	def get_SNPlist(self):
		return self.SNPlist

	def set_nlist(self, nlist):
		self.nlist = nlist

	def get_nlist(self):
		return self.nlist

	def set_genelist(self, genelist):
		self.genelist = genelist

	def get_genelist(self):
		return self.genelist

	def set_desertlist(self, desertlist):
		self.desertlist = desertlist

	def get_desertlist(self):
		return self.desertlist

	def set_desertgenes(self, desertgenes):
		self.desertgenes = desertgenes

	def get_desertgenes(self):
		return self.desertgenes

	def set_coveragelist(self, coveragelist):
		self.coveragelist = self.find_regions(sorted(coveragelist))

	def get_coveragelist(self):
		return self.coveragelist

	def set_contiglist(self, contiglist):
		self.contiglist = contiglist

	def get_contiglist(self):
		return self.contiglist
	
	def get_desertplot(self):
		return self.desertplot

	
	def get_average_SNP_spacing(self):
		i = 0 ; spacing = 0.0
		sortedSNPlist = sorted(self.SNPlist)
		spacing += sortedSNPlist[0]
		for i in range(len(sortedSNPlist)):
			if i+1 < len(sortedSNPlist):
				spacing += sortedSNPlist[i+1]-sortedSNPlist[i]
				i+=1
		spacing += self.size - sortedSNPlist[i]
		return spacing/len(self.SNPlist)
		

	def make_nlist(self, sequence):
		poslist = []; position=0
		splitlist = sequence.split('n')
		for i in range(len(splitlist)):
			if i != 0:
				if len(splitlist[i]) == 0:
					position += 1
				else:
					position += len(splitlist[i])+1
			else:
				if len(splitlist[i]) == 0:
					position += 1
				else:
					position += len(splitlist[i])
			poslist.append(position)
		poslist.pop(-1)
		if self.size != 0:
			self.nfrequency = len(poslist)/float(self.size)
		self.set_nlist(self.find_regions(sorted(poslist)))

	def find_regions(self, poslist):
		stretchlist = []; start = 0; stop = 0; stretch = False
		for i in range(len(poslist)):
			if stretch == False:
				start = poslist[i]; stop = poslist[i]
				if i+1 < len(poslist):
					if poslist[i+1] == stop+1:
						stretch = True
					else:
						stretchlist.append([start,stop])
				else:
					stretchlist.append([start,stop])
			else:
				if i+1 < len(poslist):
					if poslist[i+1] != poslist[i]+1:
						stop = poslist[i]
						stretchlist.append([start,stop])
						stretch = False
				else:
					stop = poslist[i]
					stretchlist.append([start, stop])
					stretch = False
		return stretchlist

	def search_for_desert_seeds(self, seedsize, seedremove, seedcutoff):
		seedlist = []; last = 0; seed = False
		if seedremove == False:
			for i in range(len(self.SNPlist)):
				if seed:
					if i+1 < len(self.SNPlist):
						if self.SNPlist[i+1]-self.SNPlist[i] < seedsize:
							stop = self.SNPlist[i]; seed = False
							if self.check_for_matches(start, stop, self.nlist, 0) != True:
								if self.check_for_matches(start, stop, self.nlist, seedcutoff) == False:
									seedlist.append([start, stop, stop-start, SNPs])
						if self.SNPlist[i+1]-self.SNPlist[i] > seedsize:
							SNPs.append(self.SNPlist[i])
					else:
						if self.check_for_matches(start, chromlength, self.nlist, 0) != True:
							if self.check_for_matches(start, stop, self.nlist, seedcutoff) == False:
								seedlist.append([start, chromlength, chromlength-start, SNPs])
				if not seed:
					if i+1 < len(self.SNPlist) and self.SNPlist[i+1]-self.SNPlist[i] > seedsize:
						start = self.SNPlist[i]; seed = True; SNPs = []
		else:
			for i in range(len(self.SNPlist)):
				if seed:
					if i+1 < len(self.SNPlist):
						if self.SNPlist[i+1]-self.SNPlist[i] < seedsize:
							stop = self.SNPlist[i]; seed = False
							if self.check_for_matches(start, stop, self.nlist, 0) != True:
								seedlist.append([start, stop, stop-start, SNPs])
						if self.SNPlist[i+1]-self.SNPlist[i] > seedsize:
							SNPs.append(self.SNPlist[i])
					else:
						if self.check_for_matches(start, chromlength, self.nlist, 0) != True:
							seedlist.append([start, chromlength, chromlength-start, SNPs])
				if not seed:
					if i+1 < len(self.SNPlist) and self.SNPlist[i+1]-self.SNPlist[i] > seedsize:
						start = self.SNPlist[i]; seed = True; SNPs = []
		return seedlist

	def check_for_matches(self, start, stop, stretchlist, cutoff):
		filtered = False; stretchsize = 0.0
		for stretch in stretchlist:
			if stretch[0] >= start and stretch[0] < stop:
				if stretch[1] <= stop:
					stretchsize += stretch[1]-stretch[0]
				elif stretch[1] > stop:
					stretchsize += stop-stretch[0]
			elif stretch[1] >= start and stretch[1] <= stop:
				stretchsize += stretch[1]-start
		if stretchsize/(stop-start) > cutoff:
			filtered = True
		return filtered

	def enlarge_deserts(self, sensitivity, seedlist):
		vcfframe = pd.DataFrame(self.SNPlist)
		for seed in seedlist:
			start = seed[0]; stop = seed[1]; SNPs = seed[3]; enlarge = True
			while enlarge == True:
				frequency = (vcfframe[0]>(start-500))&(vcfframe[0]<(start+500))
				if frequency.any() == True and frequency.value_counts()[True] > sensitivity:
					enlarge = False
				else:
					SNPs.append(start)
					if self.SNPlist.index(start)-1 >= 0:
						start = self.SNPlist[self.SNPlist.index(start)-1]
					else:
						start = 0; enlarge = False
			enlarge = True
			while enlarge == True:
				frequency = (vcfframe[0]>(stop-500))&(vcfframe[0]<(stop+500))
				if frequency.any() == True and frequency.value_counts()[True] > sensitivity:
					enlarge = False
				else:
					SNPs.append(stop)
					if stop != self.size and self.SNPlist.index(stop) < len(self.SNPlist):
						if self.SNPlist.index(stop)+1 < len(self.SNPlist):
							stop = self.SNPlist[self.SNPlist.index(stop)+1]
						else:
							stop = self.size; enlarge = False
					else:
						enlarge = False
			self.desertlist.append([start, stop, stop-start, SNPs])

	def intersect_deserts(self):
		for i in range(len(self.desertlist)):
			if i+1 < len(self.desertlist):
				if self.desertlist[i][1] > self.desertlist[i+1][0]:
					self.desertlist[i][3].append(self.desertlist[i][1])
					self.desertlist[i][3].append(self.desertlist[i+1][0])
					self.desertlist[i] = [self.desertlist[i][0], self.desertlist[i+1][1], self.desertlist[i+1][1]-self.desertlist[i][0], self.desertlist[i][3]+self.desertlist[i+1][3]]
					self.desertlist.pop(i+1)
				elif self.desertlist[i][1] == self.desertlist[i+1][0]:
					self.desertlist[i][3].append(self.desertlist[i][1])
					self.desertlist[i] = [self.desertlist[i][0], self.desertlist[i+1][1], self.desertlist[i+1][1]-self.desertlist[i][0], self.desertlist[i][3]+self.desertlist[i+1][3]]
					self.desertlist.pop(i+1)

	def post_process_deserts(self, minsize, desertremove, covcutoff, desertcutoff, desertfitting):
		changelist = []
		if desertfitting == True:
			self.fit_deserts_to_contigs()
		if desertremove == True and len(self.nlist) > 1:
			for desert in self.desertlist:
				if desert[2] > minsize:
					if len(self.coveragelist) > 1:
						if self.check_for_matches(desert[0], desert[1], self.coveragelist, covcutoff) == False:
							if self.check_for_matches(desert[0], desert[1], self.nlist, desertcutoff) == False:
								changelist.append(desert)
					else:
						if self.check_for_matches(desert[0], desert[1], self.nlist, desertcutoff) == False:
							changelist.append(desert)
						
		else:
			for desert in self.desertlist:
				if desert[2] > minsize:
					if len(self.coveragelist) > 1:
						if self.check_for_matches(desert[0], desert[1], self.coveragelist, covcutoff) == False:
							changelist.append(desert)
					else:
						changelist.append(desert)
		self.set_desertlist(changelist)

	def fit_deserts_to_contigs(self,):
		changelist = []
		for desert in self.desertlist:
			for contig in self.contiglist:
				if desert[0] > contig[1] and desert[1] < contig[1]:
					changelist.append(desert)
				elif desert[0] < contig[1] and desert[1] < contig[1] and desert[1] > contig[0]:
					changelist.append([contig[0],desert[1],desert[1]-contig[0],self.SNP_cutter(contig[0],desert[1],desert[3])])
				elif desert[0] > contig[0] and desert[0] < contig[1] and desert[1] > contig[1]:
					changelist.append([desert[0],contig[1],contig[1]-desert[0],self.SNP_cutter(desert[0],contig[1],desert[3])])
				elif desert[0] < contig[0] and desert[1] > contig[1]:
					changelist.append([contig[0],contig[1],contig[1]-contig[0],self.SNP_cutter(contig[0],contig[1],desert[3])])
		self.desertlist = changelist
		return None

	def SNP_cutter(self, start, stop, SNPs):
		SNPlist = []
		for SNP in SNPs:
			if SNP > start and SNP < stop:
				SNPlist.append(SNP)
		return SNPlist
		
	def find_desert_genes(self,):
		for desert in self.desertlist:
			for gene in self.genelist:
				if gene[0]>desert[0] and gene[0]<desert[1] and gene[1]<desert[1]:
					gene.append('.')
					if desert[3] > 0:
						gene = self.check_gene_for_SNPs(desert[3], gene)
					self.desertgenes.append(gene)

	def check_gene_for_SNPs(self, SNPs, gene):
		SNPlist = []; hit = False
		for SNP in SNPs:
			if SNP>gene[0] and SNP<gene[1]:
				SNPlist.append(SNP)
				gene.append(SNPlist)
		return gene
			
	def timelines(self, y, start, stop, color):
		plt.hlines(y, start, stop, color, lw=6)
		#plt.vlines(start, y+0.03, y-0.03, color, lw=2)
		#plt.vlines(stop, y+0.03, y-0.03, color, lw=2)
	
	def make_desertplot(self, plotpath):
		#fig = plt.figure(figsize=(100,5))
		for desert in self.desertlist:
			self.timelines(1, desert[0], desert[1], 'b')
		plt.yticks([1], [self.name])
		plt.ylim(0,3)
		#plt.xticks(np.arange(0, self.size, 100000))
		plt.xlim(0, self.size)
		plt.xlabel('bp')
		plt.savefig(plotpath+'_'+self.name+'.png', dpi=150)
		return None

	def get_SNP_spacing(self):
		Spacelist = []; i = 0
		Spacelist.append(str(self.SNPlist[0]))
		for i in range(len(self.SNPlist)):
			if i+1 < len(self.SNPlist):
				Spacelist.append(str(self.SNPlist[i+1]-self.SNPlist[i]))
				i+=1
		Spacelist.append(str(self.size - self.SNPlist[i]))
		return Spacelist		

	def search_for_SNP_deserts(self, sensitivity, seedsize, minsize, desertremove, desertcutoff,seedremove, seedcutoff, covcutoff, desertfitting):
		seedlist = self.search_for_desert_seeds(seedsize, seedremove, seedcutoff)
		if sensitivity > 0:
			print 'Enlarging desert seeds on '+self.get_name()+' with frequency '+str(sensitivity)
			self.enlarge_deserts(sensitivity, seedlist)
		else:
			self.set_desertlist(seedlist)
		self.intersect_deserts()
		self.post_process_deserts(minsize, desertremove, covcutoff, desertcutoff, desertfitting)

	def process_data(self, sensitivity, seedsize, minsize, desertremove, desertcutoff, seedremove, seedcutoff, covcutoff, desertfitting):
		self.search_for_SNP_deserts(sensitivity, seedsize, minsize, desertremove, desertcutoff, seedremove, seedcutoff, covcutoff, desertfitting)
		self.find_desert_genes()
