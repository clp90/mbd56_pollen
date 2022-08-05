#!/usr/bin/env python3

''' 
-------------------------
Version history:
	v.1.0: initial build	- 01/11/2021
-------------------------

'''

import sys, os, argparse, re, random, pysam, pyBigWig, math, numpy as np
from Bio.Seq import Seq
from Bio import pairwise2

if len(sys.argv) == 1:
	print("-------------------------")
	print("mikado_refine v1.0		by Colette L. Picard, 08/04/2022")
	print("-------------------------")
	print("""Usage: mikado_refine.py [options] outprefix

This is a simple script which requires most of the same inputs as the mikado 'pick' function, plus the output of 'pick' itself,
but does some extra refining esp. for detecting ncRNAs. Refining is based on total expression and on avg. 
expression along the gene from a .bw track, with the assumption that a gene's coverage should be roughly
uniform (only egregious examples are separated based on this, and it's a tunable param). 

Assumes that you have pooled all your data into a single BAM file, and follows mostly the standard Mikado
workflow, but with a few changes at the end.

Here's an example workflow:
fasta="/path/to/your/genome/fasta.fa"
fastafai="/path/to/your/genome/fasta.fa.fai"
inbam="/path/to/your/merged/bam/file.bam"

** Assumes you've already run all transcript assembly programs; we used CLASS2, cufflinks, StringTie and Trinity

# (1) Get junctions from BAM file using portcullis
portcullis full --orientation FR --strandedness firststrand --max_length 2000 -o portcullis $fasta $inbam

# (1.5) Make mikado config file mikado_config.txt:
CLASS2_transcripts.gtf	cl	True		False
cufflinks_output/transcripts.gtf	cf	True		False
stringtie_transcripts.gtf	st	True		False
trinity.gtf	tr	True		False
starting_annotations.gtf	an	True	True

# (2) initial prep and merging of annots from all sources (mikado.yaml can be found on mikado github)
mikado configure --strand-specific --reference $fasta --junctions portcullis_filtered.pass.junctions.bed --list mikado_config.txt mikado.yaml
mikado prepare --configuration mikado.yaml -o mikado.gtf -of mikado.fa

# (2.5) ORF prediction w/ Prodigal and BLAST of proteins
prodigal -i mikado.fa -g 1 -o mikado.orfs.gff3 -f gff
BLASTDB=/location/of/blast/db
cmd="blastx -max_target_seqs 5 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop\" -num_threads 8 -query mikado.fa -db blast/uniprot_sprot_plants -out mikado_prepared.blast.tsv"

# (3) finish mikado run using appropriate scoring file (here plants, plant.yaml on mikado github)
mikado serialise --transcripts mikado.fa --configuration mikado.yaml --orfs mikado.orfs.gff3 --xml mikado_prepared.blast.tsv --genome_fai $fastafai
mikado pick --configuration mikado.yaml --scoring-file plant.yaml --no-purge --output-dir mikado_pick_v1 mikado.gtf

# (3.5) Drop superloci
awk -F$'\t' '{OFS=FS} $3 != "superlocus"' mikado_pick_orig/mikado.loci.gff3 > mikado_pick_orig/mikado.loci.filt.gff3

# (4) Get bamCoverage tracks from original input BAM file (using deeptools bamCoverage):
samtools view -b -f 128 -F 16 $inbam > fwd1.bam
samtools view -b -f 80 $inbam > fwd2.bam
samtools view -b -f 144 $inbam > rev1.bam
samtools view -b -f 64 -F 16 $inbam > rev2.bam
samtools merge -f fwd.bam fwd1.bam fwd2.bam; samtools index fwd.bam
samtools merge -f rev.bam rev1.bam rev2.bam; samtools index rev.bam
bamCoverage -b fwd.bam -o fwd.bw --normalizeUsing CPM --binSize 10
bamCoverage -b rev.bam -o rev.bw --normalizeUsing CPM --binSize 10

# Once steps above complete, can run this script as a final refinement of the mikado:

mikado_refine.py mikado_pick_refined --pick mikado_pick_orig/mikado.loci.gff3 --raw orig_annots.gtf --genome $fasta --plustrack plustrack fwd.bw --minustrack rev.bw --junc portcullis_all.junctions.bed

# Note - all inputs must be sorted by start of transcript, with all entries for that transcript following
# in order (GTF/GFF3 files) - this is default output from the commands above. 
# Note 2 - mikado prepare config file should be set to put 'an' as prefix for transcripts from original annotation (as in example mikado_config.txt above), or these
# will not be properly identified (code will still run, but original annotations won't be given any additional weight)


""")
	print("-------------------------")
	sys.exit(0)

# read in arguments
parser = argparse.ArgumentParser()
parser.add_argument('outprefix', help = 'Prefix for output files')

# required inputs
parser.add_argument('--pick', default = None, help = 'Output of mikado pick (a GFF3 file)')
parser.add_argument('--raw', default = None, help = 'Output of mikado prepare (a GTF file)')
parser.add_argument('--junc', default = None, help = 'Portcullis junction BED file, output of 2nd step (*/2-junc/portcullis_all.junctions.bed)')
parser.add_argument('--genome', default = None, help = 'Genome in fasta format')
parser.add_argument('--plustrack', default = None, help = 'BigWig track of plus strand expression levels')
parser.add_argument('--minustrack', default = None, help = 'BigWig track of minus strand expression levels')

# parameters used for transcript selection
parser.add_argument('--maxlen', default = 30000, type=int, help = 'Max total length of transcript (including introns)')
parser.add_argument('--min_exon_ratio', default = 0.1, type=float, help = 'Minimum fraction of mRNA corresponding to exons (= tot len exons / tot len unspliced)')
parser.add_argument('--min_denovo_expr', default = 0.011, type=float, help = 'Minimum average expression in bigwig track normalized by CPM (average across all exons) to retain novel transcripts')
parser.add_argument('--weight_cdslen', default = 2.5, type=float, help = 'Weight given to length of CDS when selecting best transcripts')
parser.add_argument('--weight_totlen', default = 0.2, type=float, help = 'Weight given to mRNA length when selecting best transcripts')
parser.add_argument('--weight_totcov', default = 1, type=float, help = 'Weight given to total expression when selecting best transcripts')
parser.add_argument('--weight_avgcov', default = 1, type=float, help = 'Weight given to average expression when selecting best transcripts')
parser.add_argument('--weight_stdevcov', default = 0.2, type=float, help = 'Weight given to variability of expression along transcript (less is better) when selecting best transcripts')
parser.add_argument('--weight_juncreads', default = 0.75, type=float, help = 'Weight given to number of supported junction reads when selecting best transcripts')
parser.add_argument('--weight_utrdist', default = 0.2, type=float, help = 'Weight given to how far away UTR length is from Arabidopsis average (default is 155bp for 5\' UTR, 242bp for 3\' UTR, see options below)')
parser.add_argument('--weight_isorig', default = 0.15, type=float, help = 'Weight given to original annotations when selecting best transcripts')
parser.add_argument('--weight_ispicked', default = 0.1, type=float, help = 'Weight given to mikado picked transcripts when selecting best transcripts')
parser.add_argument('--weight_outside_internal', default = 0.75, type=float, help = 'Weight given to expression within gene bounds not accounted for by transcripts (less is better) when selecting best transcripts')
parser.add_argument('--weight_outside', default = 1.25, type=float, help = 'Weight given to expression outside gene bounds not accounted for by transcripts (less is better) when selecting best transcripts')

# other parameters
parser.add_argument('--min_spliced_len', default = 50, type=int, help = 'Minimum spliced length of mRNA; all transcripts below this length are ignored')
parser.add_argument('--min_prot_len', default = 50, type=int, help = 'Minimum length of protein for transcript to be considered coding')
parser.add_argument('--min_orf_frac', default = 0.35, type=float, help = 'Minimum fraction of mRNA belonging to ORF for transcript to be considered coding')
parser.add_argument('--avg_5pUTR_length', default = 155, type=int, help = 'Average length of 5\' UTR in genome')
parser.add_argument('--avg_3pUTR_length', default = 242, type=float, help = 'Average length of 3\' UTR in genome')
parser.add_argument('--fuzziness', default = 50, type=int, help = 'Allows transcripts overlapping by this many bp to still be considered nonoverlapping (edges can be poorly defined)')
parser.add_argument('--printstats', default = "", help = 'Print stats used to pick isoforms for superlocus including this transcript')

# runtime optimizing parameters (usually don't mess with these)
parser.add_argument('--max_recursion', default = 1000, type=int, help = 'Max number of gene combos identified by recursion allowed')


args = parser.parse_args()


#-------------------------------------------------------------
# HELPER FUNCTIONS
#-------------------------------------------------------------

# little function that reports an internal error which is likely due to a bug rather
# than user error, then exits
def internalError(errmsg):
	print("Internal error:",errmsg)
	print("This error is likely due to a bug in the script. Please contact script author for fix.")
	sys.exit(1)

# quick function that returns whether or not the two provides lists of tuples contain the
# same elements; input lists need not be in same order
def compareCoordList(list1,list2):
	list1sort = sorted(list1, key = lambda x: x[0])
	list2sort = sorted(list2, key = lambda x: x[0])
	if len(list1) != len(list2):
		return False
	elif list1sort != list2sort:
		return False
	else:
		return True

# quick function that returns whether or not two ranges provided as (left, right) tuples (1-based inclusive)
# overlap or not, accounting for allowed fuzziness (if overlap is < fuzziness, then they -don't- overlap)
def rangeOverlap(range1,range2):
	if (range1[0] <= range2[0] and range1[1] >= range2[0]) or (range2[0] <= range1[0] and range2[1] >= range1[0]):
		return True
	else:
		return False

# recursive function to get all valid nonoverlapping intervals, from
# StackOverflow user MBo posted Nov. 6 2018, see https://stackoverflow.com/questions/53176104/find-all-combinations-with-non-overlapped-regions
# allows up to 40 bp overlap for fuzzy overlap
def getValidIntervalsFuzzy(l, idx, right, ll, idxcombos, cont, fuzz = args.fuzziness):
	if cont == True:
		# reached base case
		if idx == len(l):
			if ll:
				idxcombos.append(ll)	
			if len(idxcombos) > args.max_recursion:
				cont = False
			return cont
		
		#find next non-overlapping interval without using l[idx]
		next = idx + 1  
		while next < len(l) and right >= l[next][0] + fuzz:
			next += 1
		cont = getValidIntervalsFuzzy(l, next, right, ll, idxcombos, cont, fuzz)

		#find next non-overlapping interval after using l[idx]
		next = idx + 1
		right = l[idx][1]
		while next < len(l) and right >= l[next][0] + fuzz:
			next += 1
		if ll == '':
			upd = ll + str(idx)
		else:
			upd = ll + ';' + str(idx)
		cont = getValidIntervalsFuzzy(l, next, right, upd, idxcombos, cont, fuzz)

		return cont



#-------------------------------------------------------------
# CLASSES
#-------------------------------------------------------------

# Transcript 
# --------------------
# represents one transcript and stores all relevant info, including
# general features (chromosome, name, source, etc.), feature coordinates 
# (locations of exons, CDS etc.), and features used for scoring (expression 
# data, original status etc.)
class Transcript:

	# -------------
	# initialize new Transcript object from name, chromosome, strand and source
	def __init__(self, tname, chromosome, strand, source):
		
		# basic info about this transcript
		self.origname = tname
		self.finalname = None
		self.chr = chromosome
		self.strand = strand
		self.source = source
		self.winner = False
		self.finished = False
		self.is_antisense = False
		self.forcedCDS = False
		
		# coordinates of gene features (exons, UTRs, CDSs etc.) all saved as lists of tuples
		self.mRNA = ()		# overall transcript bounds
		self.exons = []
		self.CDS = []
		self.UTR5p = []
		self.UTR3p = []
		
		# other transcript features used for scoring
		self.orig_annot = tname.startswith("an_")	# annotation was present in original annotations
		self.picked = False							# whether or not this transcript was picked by Mikado
		self.totCDS = None							# total length of CDS (0 if noncoding)
		self.totlen = None							# total spliced length of mRNA
		self.coding = False							# contains valid coding sequence
		self.tot_supp_juncs = 0						# total junction reads overlapping this transcript (regardless of whether consistent w/ exons)
		self.frac_supp_juncs = 0					# fraction of total junction reads consistent w/ this transcript's exons
		self.var_all = 100							# avg. coverage variability along entire transcript (omits 5% of each end); 100 by default as this should be higher than possible values
		self.avg_cov = 0							# avg. coverage overall
		self.tot_cov = 0							# total expression over this locus
		self.UTR_dist = -1
		self.protein = ""							# protein sequence of this transcript
		
	# -------------
	# basic function to add new info from a line from GTF file; only intended for mRNA, transcript, exon and CDS lines
	def addLine(self, type, istart, iend):
		if iend > istart:			# apparently it does need to be said, that an interval must contain at least one base...
			if type == "mRNA" or type == "transcript":
				if len(self.mRNA) == 0:
					self.mRNA = (istart,iend)
				else:
					if self.mRNA[0] != istart or self.mRNA[1] != iend:
						print("Error: transcript",self.origname,"has two lines with feature type (3rd column) == \"mRNA\" or \"transcript\" that have different coordinates")
						sys.exit(1)
			elif type == "exon":
				# append to list and resort
				self.exons.append((istart,iend))
			elif type == "CDS":
				self.CDS.append((istart,iend))
			else:
				# unrecognized type
				errstr = "tried to use addLine() with type ="+type
				internalError(errstr)
		else:
			print("Warning: GTF line has invalid interval:",self.chr+":"+str(istart)+"-"+str(iend))	

	# -------------
	# once all exons and CDS entries are added, this adds in the 5' and 3' UTRs (if protein-coding)
	def addUTRs(self):
		if self.exons == []:
			errstr = "tried to run addUTRs() on transcript without any exons, transcript name ="+self.origname
			internalError(errstr)
		if len(self.CDS) > 0:

			leftUTR = []; rightUTR = []

			# get leftmost CDS endpoint
			leftmostCDS = self.CDS[0][0]
			rightmostCDS = self.CDS[-1][1]
			
			# start at first exon and go thru all exons with right endpoint before 
			# or equal to leftmostCDS to get leftUTR
			i=0
			while i < len(self.exons) and self.exons[i][0] < leftmostCDS:
#				print("exon",self.exons[i],"still less than CDS")
				if self.exons[i][1] >= leftmostCDS:
					leftUTR.append((self.exons[i][0], leftmostCDS-1))
				else:
					leftUTR.append((self.exons[i][0], self.exons[i][1]))
				i+=1
			
			# start at last exon and go backwards to get the rightUTR
			i=len(self.exons)-1
			while i >= 0 and self.exons[i][1] > rightmostCDS:
#				print("exon",self.exons[i],"still more than CDS")
				if self.exons[i][0] <= rightmostCDS:
					rightUTR.append((rightmostCDS+1, self.exons[i][1]))
				else:
					rightUTR.append((self.exons[i][0], self.exons[i][1]))
				i-=1
			
			# update UTRs
			if self.strand == "+":
				self.UTR5p = leftUTR; self.UTR3p = rightUTR
			else:
				self.UTR5p = rightUTR; self.UTR3p = leftUTR
				
	# -------------
	# once all exons and CDS entries are added, this adds the overall bounds
	def updatemRNA(self):		
		mRNAcoord = (min(i for (i,j) in self.exons), max(j for (i,j) in self.exons))
		# if mRNA bounds were already provided in GTF file, check they're the same
		if len(self.mRNA) != 0:
			if self.mRNA[0] != mRNAcoord[0] or self.mRNA[1] != mRNAcoord[1]:
				print("Warning: mRNA boundaries identified for transcript",self.origname,"manually don't match the one in --raw GTF file:")
				print(self.mRNA)
				print("Using identified bounds instead:",mRNAcoord)
		self.mRNA = mRNAcoord			
		
	
	# -------------
	# adds junction information, both total overlapping junctions and supported junctions
	# to calculate fraction supported junctions
	def addJuncs(self):
		totoverlapping = 0
		totsupported = 0
		
		# list all possible junctions for this transcript (assumes exons sorted)
		possiblejuncs = []
		for i in range(1,len(self.exons)):
			possiblejuncs.append((self.exons[i-1][1],self.exons[i][0]))
		
		# get all junctions that overlap this transcript
	#	print("Getting all junctions in interval",self.mRNA,"on",self.chr,self.strand,"strand")
		for kkey, terms in juncs[self.chr+self.strand].items():
			if kkey > self.mRNA[0] and kkey < self.mRNA[1]:
				for tt in terms:
	#				print("Junction",kkey,"-",tt[0],"with depth",tt[1],"overlaps interval")
					totoverlapping += tt[1]
					
					# is this junction consistent with this transcript?
					if (kkey,tt[0]) in possiblejuncs:
	#					print("Junction is consistent with transcript, adding",tt[1],"to supported count")
						totsupported += tt[1]
			elif kkey > self.mRNA[1]:
				# out of bounds, stop loop
				break

		self.tot_supp_juncs = totsupported
		if totoverlapping > 0:
			self.frac_supp_juncs = totsupported/totoverlapping
		
	# -------------
	# checks that the CDS provided produces a valid protein starting with M and ending with */STOP
	def validateCDS(self, force = False):
		cdsseq = Seq(''.join([ genomeseq.fetch(self.chr, ex[0]-1, ex[1]) for ex in self.CDS ]))
		if self.strand == '-':
			cdsseq = cdsseq.reverse_complement()
		if len(cdsseq) % 3 != 0:
	#		print("Not a multiple of 3")
			return False,None
		prot = Seq.translate(cdsseq, table='Standard', stop_symbol='*', to_stop=False, cds=False)
		if prot[0] == "M" and prot[-1] == "*" and "*" not in prot[:-1]:
			# ORF is valid, is it long enough?
	#		print("ORF is valid")
			if self.orig_annot or force == True or (len(prot) > args.min_prot_len and ((len(prot) * 3) / (self.mRNA[1]-self.mRNA[0]+1)) > args.min_orf_frac):
	#			print("ORF passes thresholds")
				return True,prot[0:-1]
			else:
	#			print("Protein stats fail:")
	#			print("Protein length:",len(prot),"(min is",args.min_prot_len,")")
	#			print("Pct of mRNA is ORF:",((len(prot) * 3) / (self.mRNA[1]-self.mRNA[0]+1)),"(min is",args.min_orf_frac,")")
				return False,None
		else:
			return False,None		
		
	# -------------
	# adds CDS if not already provided, else validates existing CDS
	def addCDS(self, force = False):
		# if CDS already provided, ensure it's correct (starts w/ M and ends w/ *); if yes stop here
	#	print("Getting CDS for",self.origname)
		# if force is true, it will always consider the CDS valid regardless of protein length or min_orf_frac
				
		predictORFs = True
		if len(self.CDS) > 0:
			res = self.validateCDS()
			if res[0] == True:
				predictORFs = False
				self.coding = True
				self.protein = res[1]
				
		# this gene has no or faulty CDS provided, re-do CDS discovery here
		if predictORFs == True:
	
			mRNAseq = Seq(''.join([ genomeseq.fetch(self.chr, ex[0]-1, ex[1]) for ex in self.exons ]))
			if self.strand == '-':
				mRNAseq = mRNAseq.reverse_complement()
		
			# get all three frames 0,1 and 2
			frames = [Seq.translate(mRNAseq[i:]+Seq('N'*(3-(len(mRNAseq[i:]) % 3))), table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)]
							
			# find longest complete orf in each frame
			longestlen = 0; longestidx = (); longestframe = -1; idx = 0
			for frame in frames:
	#			print("Looking for ORFs in: ",frame)
	#			print(re.findall(r'M[^*]+\*', str(frame)))
	#			print([len(ll) for ll in re.findall(r'M[^*]+\*', str(frame))])
				finds = re.findall(r'M[^*]+\*', str(frame))
				if len(finds) > 0:
					longest = max(re.findall(r'M[^*]+\*', str(frame)), key = len)
					if len(longest) > longestlen:	
						startidx = frame.find(longest) + 1
						endidx = startidx + len(longest) - 1
				
						longestlen = len(longest)
						longestidx = (startidx, endidx)
						longestframe = idx
				idx += 1
				
	#		if self.origname == "cf_AT1G11530.1":
	#			print("Longest ORF was:",longestlen,"aa in frame:",longestframe)
	#			print("Index:",longestidx)			
	#			print("Protein sequence:",frames[longestframe][longestidx[0]-1:longestidx[1]])
	#			print(longestlen)
	#			print((longestlen*3) / len(mRNAseq))
			
	#			sys.exit(1)
														
			# keep this ORF?
			CDSlist = []
			if force == True or (longestlen > args.min_prot_len and (longestlen*3) / len(mRNAseq) > args.min_orf_frac):
			
				# save protein sequence
				self.protein = frames[longestframe][longestidx[0]-1:longestidx[1]-1]
			
				# translate ORF coordinates back to mRNA coordinates
				orfstart = 3*(longestidx[0]-1) + longestframe + 1
				orfend = 3*(longestidx[1]) + longestframe
				
	#			print("mRNA index:",(orfstart,orfend))
	#			print("mRNA:",mRNAseq[orfstart-1:orfend])
	#			print(self.protein)
										
				# use exon coordinates to translate overall ORF boundaries into CDS entries
				tot=0; startex = 0; endex = 0		# tot keeps running sum of bps passed in exons, startex flags if starting exon has been found, endex flags if ending exon has been found
				if self.strand == "-":
					for ex in reversed(self.exons):
	#					print("Checking overlap w/ exon:",ex)
						if tot + (ex[1]-ex[0]+1) < orfstart:
	#						print("CDS not in exon",ex,"cumulative:",tot,"to",tot+ (ex[1]-ex[0]+1))
							pass 		# doesn't do anything, this if is just useful for debugging
						elif tot + (ex[1]-ex[0]+1) >= orfstart and tot + (ex[1]-ex[0]+1) < orfend:
							if startex == 0:		# starts here but doesn't end here
								startex = 1
								CDSlist.append((ex[0],ex[1]-(orfstart-tot-1)))
	#							print("CDS starts in exon",ex,"cumulative:",tot,"to",tot+ (ex[1]-ex[0]+1))
							elif startex > 0 and endex == 0:
								CDSlist.append(ex)
	#							print("CDS fully contains exon",ex,"cumulative:",tot,"to",tot+ (ex[1]-ex[0]+1))
						else:
							if startex == 0 and endex == 0:		# starts and ends in same exon
	#							print("CDS begins and ends in exon",ex,":",tot,"to",tot+ (ex[1]-ex[0]+1))
								startex = 1; endex = 1
								CDSlist.append((ex[1]-(orfend-tot-1),ex[1]-(orfstart-tot-1)))
							elif endex == 0:		# ends here
	#							print("CDS ends in this exon:",tot,"to",tot+ (ex[1]-ex[0]+1))
								endex = 1
								CDSlist.append((ex[1]-(orfend-tot-1),ex[1]))
							else:
	#							print("CDS not present in exon",ex,":",tot,"to",tot+ (ex[1]-ex[0]+1))
								pass
						tot += (ex[1]-ex[0]+1)
					CDSlist = CDSlist[::-1]
				else:
					for ex in self.exons:
	#					print("Checking overlap w/ exon:",ex)
						if tot + (ex[1]-ex[0]+1) < orfstart:
	#						print("CDS not in exon",ex,"cumulative:",tot,"to",tot+ (ex[1]-ex[0]+1))
							pass 		# doesn't do anything, this if is just useful for debugging
						elif tot + (ex[1]-ex[0]+1) >= orfstart and tot + (ex[1]-ex[0]+1) < orfend:
							if startex == 0:		# starts here but doesn't end here
								startex = 1
								CDSlist.append((ex[0]+(orfstart-tot)-1,ex[1]))
	#							print("CDS starts in exon",ex,"cumulative:",tot,"to",tot+ (ex[1]-ex[0]+1))
							elif startex > 0 and endex == 0:
								CDSlist.append(ex)
	#							print("CDS fully contains exon",ex,"cumulative:",tot,"to",tot+ (ex[1]-ex[0]+1))
						else:
							if startex == 0 and endex == 0:		# starts and ends in same exon
	#							print("CDS begins and ends in exon",ex,":",tot,"to",tot+ (ex[1]-ex[0]+1))
								startex = 1; endex = 1
								CDSlist.append((ex[0]+(orfstart-tot)-1,ex[0]+(orfend-tot)-1))
							elif endex == 0:		# ends here
	#							print("CDS ends in exon",ex,":",tot,"to",tot+ (ex[1]-ex[0]+1))
								endex = 1
								CDSlist.append((ex[0],ex[0]+(orfend-tot)-1))
							else:
	#							print("CDS not present in exon",ex,":",tot,"to",tot+ (ex[1]-ex[0]+1))
								pass
						tot += (ex[1]-ex[0]+1)
				
	#			print("List of exons:",self.exons)					
	#			print("List of CDSs:",CDSlist)
				
				self.CDS = CDSlist
				self.coding = True
				if force ==  True:
					self.forcedCDS = True
			
			else:
				# no good ORFs, gene is noncoding
				self.CDS = []
				self.coding = False

		
	# -------------
	# once all exons, etc. added, finish up everything by sorting, updating fields etc.
	def finishFullEntry(self):
	#	print("Finishing",self.origname)
		self.finished = True
	
		# clean up exons, get mRNA boundaries
		self.exons = sorted(list(set(self.exons)), key = lambda x: x[0])		# drop redundancies if any
		self.updatemRNA()
		self.totlen = sum([ex[1]-ex[0]+1 for ex in self.exons])
		
		# don't bother with the rest if this isn't a valid transcript
		if self.isValid():
		
			# add junction support
			self.addJuncs()
				
			# add rough measure of coverage and coverage variability across exons from bigwig files, if provided
			# for std dev, cut 5% of spliced length off edges to avoid edge differences; at least 10bp
			leftedge = max(10,round(self.totlen / 100 * 5))
			rightedge = self.totlen - leftedge
			tot = 0; centralcoord = []
			
			for ex in self.exons:
				cur = ex[1] - ex[0] + 1
	#			print("Exon is",ex,"and current total is:",tot)
				if tot + cur >= leftedge and tot + cur < leftedge + (ex[1]-ex[0]+1):	
	#				print("Left edge is in this exon")
					centralcoord.append((ex[0]+(leftedge-tot),ex[1]))
				elif tot + cur >= leftedge and tot + cur < rightedge:
	#				print("Fully contains this exon")
					centralcoord.append(ex)
				elif tot + cur >= rightedge:
	#				print("Right edge in this exon")
					centralcoord.append((ex[0],ex[0]+(rightedge-tot)))
					break
				tot += cur

			exons_touse = [ex for ex in self.exons if ex[1]-ex[0]+1 >= 20]
			centralcoord_touse = [ex for ex in centralcoord if ex[1]-ex[0]+1 >= 20]
			
			if len(exons_touse) > 0:
				avgcov = [None]; totcov = [None]
				if self.strand == "+" and not args.plustrack is None:
					avgcov = [plustrack.stats(self.chr, ex[0], ex[1])[0] for ex in exons_touse]		# note - exon length must be > binsize of bigwig or will fail; I just set this b/c I usually use binsize 1 or 10
					totcov = [plustrack.stats(self.chr, ex[0], ex[1], type = "sum", exact=True)[0] for ex in exons_touse]		# note - exon length must be > binsize of bigwig or will fail; I just set this b/c I usually use binsize 1 or 10
				elif not args.minustrack is None:
					avgcov = [minustrack.stats(self.chr, ex[0], ex[1])[0] for ex in exons_touse]
					totcov = [minustrack.stats(self.chr, ex[0], ex[1], type = "sum", exact=True)[0] for ex in exons_touse]		# note - exon length must be > binsize of bigwig or will fail; I just set this b/c I usually use binsize 1 or 10
								
				# convert result to weighted average
				if not (len(avgcov) == 1 and avgcov[0] is None) and len(avgcov) > 0:				
					# get weighted avg coverage (weight by total exon length)
					self.avg_cov = 0; tot = 0
	#				print("Name:",self.origname)
	#				print("exons:",exons_touse,"; avg cov:",avgcov)
					for i in range(len(avgcov)):
						self.avg_cov += avgcov[i]*(exons_touse[i][1] - exons_touse[i][0] + 1)
						tot += exons_touse[i][1] - exons_touse[i][0] + 1
					self.avg_cov = self.avg_cov / tot
	#				print("Overall avg cov:",self.avg_cov)
				
				if not (len(totcov) == 1 and totcov[0] is None) and len(totcov) > 0:
					totcov = [x if x is not None else 0 for x in totcov]	
	#				print("Totcov:",totcov)			
					self.tot_cov = sum(totcov)
	#				print("Overall tot cov:",self.tot_cov)
									
			if len(centralcoord_touse) > 0:
				stdcov = [None]
				if self.strand == "+" and not args.plustrack is None:
					stdcov = [plustrack.stats(self.chr, ex[0], ex[1], type="std")[0] for ex in centralcoord_touse]
				elif not args.minustrack is None:
					stdcov = [minustrack.stats(self.chr, ex[0], ex[1], type="std")[0] for ex in centralcoord_touse]
				if not (len(stdcov) == 1 and stdcov[0] is None) and len(stdcov) > 0:
					# get measure of range of variability
					self.var_all = max(stdcov) - min(stdcov)
								
			# add CDS info if not already present and correct, then update UTRs
			self.addCDS()
			self.updateCDS()
						
			# add whether or not was picked by mikado
			if self.chr+self.strand in mikado_picked:
				if self.mRNA in mikado_picked[self.chr+self.strand]:
					for mm in mikado_picked[self.chr+self.strand][self.mRNA]:
						if self.compareExons(mm) == True:
							self.picked = True
	
	
	# -------------
	# updates everything after a new CDS is added
	def updateCDS(self):
		self.CDS = sorted(list(set(self.CDS)), key = lambda x: x[0])
		self.addUTRs()
		self.totCDS = sum([cds[1]-cds[0]+1 for cds in self.CDS])
		
		# get average distance of UTRs from Arabidopsis average (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5828884/ Fig. 1B)
		if len(self.UTR5p) > 0:
			UTR5p_len = sum([ex[1]-ex[0]+1 for ex in self.UTR5p])
			dist1 = abs(UTR5p_len-args.avg_5pUTR_length)
		if len(self.UTR3p) > 0:
			UTR3p_len = sum([ex[1]-ex[0]+1 for ex in self.UTR3p])
			dist2 = abs(UTR3p_len-args.avg_3pUTR_length)

		if len(self.UTR5p) > 0 and len(self.UTR3p) > 0:
			self.UTR_dist = (dist1+dist2)/2
		elif len(self.UTR5p) > 0:
			self.UTR_dist = (dist1+args.avg_3pUTR_length)/2
		elif len(self.UTR3p) > 0:
			self.UTR_dist = (dist2+args.avg_5pUTR_length)/2
		else:
			if self.coding:
				self.UTR_dist = (args.avg_3pUTR_length+args.avg_5pUTR_length)/2
			# if noncoding, it stays as -1
			
		
	# -------------
	# minimal version of finishFullEntry() just to sort exons and update gene bounds
	def finishPartialEntry(self):
		# clean up exons, get mRNA boundaries
		self.exons = sorted(list(set(self.exons)), key = lambda x: x[0])		# drop redundancies if any
		self.totlen = sum([ex[1]-ex[0]+1 for ex in self.exons])
		self.updatemRNA()
		
	# -------------
	# compares two transcripts across all exons and CDS intervals, and outputs similarity score
	# based on result in both numeric and string format
	def compareTranscripts(self, tr2):
		# returns a numeric code and a human-interpretable string comparing
		# this transcript object to another one coordinate-wise
		
		# If transcript is coding:
		# 0 = perfect match on all feature coordinates
		# 1 = perfect match on CDS, one UTR different
		# 2 = perfect match on CDS, both UTRs different
		# 3 = altered CDS, but greater than 99% identity
		# 4 = altered CDS, greater than 95% identity
		# 5 = altered CDS, greater than 80% identity
		# 6 = altered CDS, greater than 50% identity
		# 7 = some overlap but < 50% identity of CDS
		# 8 = no overlap
		
		# If noncoding:
		# 0 = perfect match on all exons
		# 1 = altered exons with very high overlap (~ >= 95%)
		# 2 = altered exons with high overlap (80-95%)
		# 3 = altered exons with medium overlap (~ 50-80%)
		# 4 = altered exons with low overlap (<50%)
		# 5 = no overlap
		
		# If one is coding and one noncoding, either:
		# 7 = some overlap but < 50% identity of CDS
		# 8 = no overlap
		
				
		# sanity check
		if self.chr != tr2.chr:
			print("Internal Error: tried to compare two transcripts not on same chromosome; something went wrong with superlocus assignment")
			print(self.isoStats)
			print(tr2.isoStats)
			print("This error is likely due to a bug in the script. Please contact script author for fix.")
			sys.exit(1)

		# compare the two transcripts
		if self.coding == True and tr2.coding == True:
			if max(0, min(self.mRNA[1], tr2.mRNA[1]) - max(self.mRNA[0], tr2.mRNA[0])) <= 0:
				return 8,"no overlap"
			else:
				exonmatch = compareCoordList(self.exons,tr2.exons)
				CDSmatch = compareCoordList(self.CDS, tr2.CDS)
				utr5match = compareCoordList(self.UTR5p, tr2.UTR5p)
				utr3match = compareCoordList(self.UTR3p, tr2.UTR3p)
				
				if exonmatch == True:
					return 0,"perfect match"
				elif CDSmatch == True and ((utr5match == True and utr3match == False) or (utr5match == False and utr3match == True)):
					return 1,"CDS match, one UTR different"
				elif CDSmatch == True and utr5match == False and utr3match == False:
					return 2,"CDS match, both UTRs different"
				else:
					# note - identity will be low if frameshift! That's by design...
		#			print(self.origname,"protein:",self.protein)
		#			print(tr2.origname,"protein:",tr2.protein)
					matched = pairwise2.align.globalxx(self.protein, tr2.protein, one_alignment_only = True)[0]
		#			print(matched)
					pct_identity = (matched[2] / matched[4]) * 100
					
					if pct_identity >= 99:
						return 3,"altered CDS with >=99% identity"
					elif pct_identity >= 95:
						return 4,"altered CDS with 95-99% identity"
					elif pct_identity >= 80:
						return 5,"altered CDS with 80-95% identity"
					elif pct_identity >= 50:
						return 6,"altered CDS with 50%-80% identity"
					else:
						return 7,"altered CDS with <50% identity"

		elif self.coding == False and tr2.coding == False:
			if max(0, min(self.mRNA[1], tr2.mRNA[1]) - max(self.mRNA[0], tr2.mRNA[0])) <= 0:
				return 5,"no overlap"
			else:
				exonmatch = compareCoordList(self.exons,tr2.exons)		
				if exonmatch == True:
					return 0,"perfect match"
				else:
					# Thanks to TemporalWolf on StackOverflow edited Jun 2018 for this very compact solution
					# https://stackoverflow.com/questions/40367461/intersection-of-two-lists-of-ranges-in-python/40368603
					overlap = [[max(first[0], second[0]), min(first[1], second[1])] for first in self.exons for second in tr2.exons if max(first[0], second[0]) <= min(first[1], second[1])]
					lenoverlap = sum([a[1]-a[0]+1 for a in overlap])
					
					# take average of percent overlap from both transcripts
					pct_overlap = (((lenoverlap / self.totlen) + (lenoverlap / tr2.totlen)) / 2) * 100
					
					if pct_overlap >= 95:
						return 1,"altered exons with very high overlap (~ >= 95%)"
					elif pct_overlap >= 80:
						return 2,"altered exons with high overlap (80-95%)"
					elif pct_overlap >= 50:
						return 3,"altered exons with medium overlap (~ 50-80%)"
					else:
						return 4,"altered exons with low overlap (<50%)"
		else:
			# if one is coding and the other is noncoding but they overlap somewhat, return 7 (CDS no match)
			# if overlap, else 8 (no overlap)
			if max(0, min(self.mRNA[1], tr2.mRNA[1]) - max(self.mRNA[0], tr2.mRNA[0])) <= 0:
				return 8,"no overlap"
			else:
				return 7,"some overlap but different coding status"
		
			
	# -------------
	# simpler version of compareTranscripts() that just outputs whether the set of exons for two transcripts is the same
	def compareExons(self, tr2):
		if self.chr != tr2.chr or self.strand != tr2.strand:
			return False
		return compareCoordList(self.exons,tr2.exons)
	
	# -------------
	# returns how much two transcripts overlap using only their boundaries, ignoring exons
	def overlapping(self,tr2):
		if self.chr == tr2.chr:
			return max(0, min(self.mRNA[1], tr2.mRNA[1]) - max(self.mRNA[0], tr2.mRNA[0]))
		else:
			return 0
			
	# -------------
	# same as overlapping(), but expects argument to be a superlocus instead of another transcript
	def overlapping2(self,superlocus):
		if self.chr == superlocus.chr:
			return max(0, min(self.mRNA[1], superlocus.boundaries[1]) - max(self.mRNA[0], superlocus.boundaries[0]))
		else:
			return 0
			
	# -------------
	# returns how many bp overlap between the exons of two transcripts
	# currently brute forced, could optimize but these lists are typically very short, probably not worth speeding up
	def overlap_by_exon(self,tr2):
		totoverlap = 0
		if self.chr == tr2.chr:
			for e1 in self.exons:
				for e2 in tr2.exons:
					totoverlap += max(0, min(e1[1], e2[1]) - max(e1[0], e2[0])) + 1
		return totoverlap	
		
	# -------------
	# same as above, but expects a list of tuples corresponding to exons, instead of a Transcript object
	def overlap_by_exon2(self,exonlist):
		totoverlap = 0
		for e1 in self.exons:
			for e2 in exonlist:
				totoverlap += max(0, min(e1[1], e2[1]) - max(e1[0], e2[0])) + 1
		return totoverlap	

	# -------------
	# outputs true/false whether this transcript is valid based on several factors
	# returns tuple of true/false and string explanation
	def isValid(self):
		self.finishPartialEntry()
		if self.totlen < args.min_spliced_len:
			return False
		exonratio = self.totlen / (self.mRNA[1]-self.mRNA[0]+1)
		if self.strand not in ['+','-']:
			return False
		elif self.mRNA[1] - self.mRNA[0] > args.maxlen:
			return False
		elif exonratio < args.min_exon_ratio:
			return False
		else:
			return True
			
	# -------------
	# write this transcript to file in GTF format
	def writeToGTF(self, outfile):
		gene_id = self.finalname		
		towrite = [self.chr,self.source,"gene",str(self.mRNA[0]),str(self.mRNA[1]),".",self.strand,".",'gene_id "'+gene_id+'"; transcript_id "'+gene_id+'"\n']
		outfile.write('\t'.join(towrite))
		towrite = [self.chr,self.source,"mRNA",str(self.mRNA[0]),str(self.mRNA[1]),".",self.strand,".",'gene_id "'+gene_id+'"; transcript_id "'+gene_id+'"\n']
		outfile.write('\t'.join(towrite))
		if self.UTR5p != None and self.UTR5p != []:
			for u in sorted(self.UTR5p, key = lambda x: x[0]):
				towrite = [self.chr,self.source,"five_prime_UTR",str(u[0]),str(u[1]),".",self.strand,".",'gene_id "'+gene_id+'"; transcript_id "'+gene_id+'"\n']
				outfile.write('\t'.join(towrite))
		
		for exon in sorted(self.exons, key = lambda x: x[0]):
			towrite = [self.chr,self.source,"exon",str(exon[0]),str(exon[1]),".",self.strand,".",'gene_id "'+gene_id+'"; transcript_id "'+gene_id+'"\n']
			outfile.write('\t'.join(towrite))
			
		for cds in sorted(self.CDS, key = lambda x: x[0]):
			towrite = [self.chr,self.source,"CDS",str(cds[0]),str(cds[1]),".",self.strand,".",'gene_id "'+gene_id+'"; transcript_id "'+gene_id+'"\n']
			outfile.write('\t'.join(towrite))
		
		if self.UTR3p != None and self.UTR3p != []:
			for u in sorted(self.UTR3p, key = lambda x: x[0]):
				towrite = [self.chr,self.source,"three_prime_UTR",str(u[0]),str(u[1]),".",self.strand,".",'gene_id "'+gene_id+'"; transcript_id "'+gene_id+'"\n']
				outfile.write('\t'.join(towrite))
				
	# -------------
	# print info for this transcript
	def printTranscript(self):
		# prints useful info about this transcript
		print("----------")
		if self.winner == True:
			print("WINNER")
		print("Original transcript name:",self.origname)
		if not self.finalname is None:
			print("Final transcript name:",self.finalname)
		print("----------")
		print(self.source,"transcript on",self.chr,self.strand,"strand")
		print("Original annotation?",self.orig_annot)
		print("Picked annotation?",self.picked)
		print("Total spliced length:",self.totlen)
		print("Coding?",self.coding)
		if self.coding:
			print(" - CDS len:",self.totCDS)
		print("Avg coverage:",self.avg_cov)
		print("Cov variability along middle 90% of mRNA:",self.var_all)
		print("Total junction reads supporting transcript:",self.tot_supp_juncs)
		if self.tot_supp_juncs != 0:
			print("Total fraction of junction reads supporting transcript:",self.frac_supp_juncs)
		else:
			print("No junction reads detected")
		print("----------")
		print("Overall boundaries:",self.mRNA)
		print("Exons:",self.exons)
		if self.coding:
			print("CDS:",self.CDS)
			print("Protein sequence:",self.protein)
		else:
			print("Noncoding")	

			
	# -------------
	# avoid collisions in naming by making sure name has not already been assigned
	def addName(self, newname):
	
		global usednames
		
		if newname in usednames:
			usednames[newname] += 1
			self.finalname = newname + '_v' + str(usednames[newname])
		else:
			usednames[newname] = 1
			self.finalname = newname
		
	
# Superlocus 
# --------------------
# represents a group of multiple transcripts overlapping at least to some degree
# on the same strand. Can include transcripts that don't overlap with each other
# but overlap with other members of the superlocus.
class Superlocus:
	
	def __init__(self, transcript):
		
		# basic info about this superlocus
		self.transcripts = [transcript]
		self.chr = transcript.chr
		self.strand = transcript.strand
		self.boundaries = (transcript.mRNA[0],transcript.mRNA[1])
	
	def addTranscript(self, transcript):
		# add new Transcript, after checking if it's a duplicate of an existing Transcript
		isdup = False
		for i in range(len(self.transcripts)):
			if self.transcripts[i].compareExons(transcript) == True:
				# this exact transcript is already there, dont' add it
				isdup = True
				# if the new transcript is an original annotation, then replace the old entry
				if transcript.orig_annot == True:
					self.transcripts[i] = transcript
				break
		if isdup is False:
			# transcript isn't already present, add it to superlocus
	#		print("Not a duplicate, adding",transcript.origname,"to superlociplus")
			self.transcripts.append(transcript)
			# extend boundaries of this superlocus
			self.boundaries = (min(self.boundaries[0],transcript.mRNA[0]),max(self.boundaries[1],transcript.mRNA[1]))	
	
	def finishSuperlocus(self):
		# for all transcripts in this superlocus, if one is an original annot with a CDS, and others
		# in the superlocus fully contain all the same CDSs as the original annot but were considered
		# noncoding due to not meeting thresholds, force add the CDSs back so that they're all considered coding
	#	print("Finishing superlocus containing:",[x.origname for x in self.transcripts])
		origannots = [x for x in self.transcripts if x.orig_annot == True and x.coding == True]				# coding original annots
		novelannots = [x for x in self.transcripts if x.orig_annot == False]								# noncoding original annots
		if len(origannots) > 0 and len(novelannots) > 0:
			for tt in novelannots:
				for oo in origannots:
					if tt.coding == False:				
	#					print("Checking if novel annot",tt.origname,"overlaps CDS of original",oo.origname)
	#					print("Novel annot exons:",tt.exons,"and original CDS:",oo.CDS)
						if tt.overlap_by_exon2(oo.CDS)/oo.totCDS == 1:
	#						print("Exons of novel transcript",tt.origname,"contain CDS of original",oo.origname)
							tt.addCDS(force = True)
							tt.updateCDS()
	#					else:
	#						print("Overlap is",tt.overlap_by_exon2(oo.CDS),"and total CDS len=",oo.totCDS)
					else:
	#					print("Checking if novel coding annot",tt.origname,"has same CDS as original",oo.origname)
	#					print("Novel annot CDS:",tt.CDS,"and original CDS:",oo.CDS)
						if compareCoordList(oo.CDS,oo.CDS) == True:
	#						print("Novel transcript",tt.origname,"has same CDS as original",oo.origname)
							tt.forcedCDS = True
						
		
	def printInfo(self):
		print("Superlocus on",self.chr,self.strand,"with boundaries ("+str(self.boundaries[0])+"-"+str(self.boundaries[1])+") and containing",len(self.transcripts),"transcripts:")
		for tt in self.transcripts:
			print(tt.origname)
				

	
	


#-------------------------------------------------------------
# MAJOR FUNCTIONS
#-------------------------------------------------------------

# globals used in major functions pickBestTranscript() and outputResults() below

# keep count of number of transcripts of each type to output in final summary
coding_stats = [0,0,0,0]				# levels, in order: identical, very similar to existing, some similarity, novel
coding_mapping = {0 : 0, 1 : 1, 2 : 1, 3 : 1, 4 : 2, 5 : 2, 6 : 3, 7 : 3, 8 : 3}
noncoding_stats = [0,0,0,0]				# same levels as coding_stats
noncoding_mapping = {0 : 0, 1 : 1, 2 : 2, 3 : 2, 4 : 3, 5 : 3}						
merged = 0
ssplit = 0

# list of all gene names that have been assigned, to avoid repeats
usednames = {}		# dict counting number of times a name has been assigned

def outputSuperlocus(superlocus):

	global coding_stats
	global noncoding_stats
	global merged
	global ssplit
	
#	print("Outputting superlocus with the following transcripts:")
#	for tt in superlocus.transcripts:
#		print(tt.origname)
	
	tooutput = []
	if len(superlocus.transcripts) == 1:
		# edge case: only one transcript here; just output it unless it's novel AND very poorly expressed AND it's noncoding
		if superlocus.transcripts[0].orig_annot == True:
			superlocus.transcripts[0].addName(superlocus.transcripts[0].origname[3:].split('.')[0])
			if superlocus.transcripts[0].coding == True:
				coding_stats[0] += 1
			else:
				noncoding_stats[0] += 1
		elif superlocus.transcripts[0].avg_cov >= args.min_denovo_expr or superlocus.transcripts[0].coding == True:
			# only keep novel transcripts with some coverage
			if superlocus.transcripts[0].coding == True:
				superlocus.transcripts[0].addName("novel_coding_"+str(coding_stats[3]))
				coding_stats[3] += 1
			else:
				superlocus.transcripts[0].addName("novel_noncoding_"+str(noncoding_stats[3]))
				noncoding_stats[3] += 1
		else:
	#		print("Novel annotation",superlocus.transcripts[0].origname,"has nearly no coverage, and was omitted")
	#		print("Coverage of the novel annot was:",superlocus.transcripts[0].avg_cov)
			superlocus.transcripts[0].winner = False
	else:
		# assign names to all winners in this superlocus
		winners = [tr for tr in superlocus.transcripts if tr.winner]
		origtranscripts = [tr for tr in superlocus.transcripts if tr.orig_annot]
		origtranscriptnames = [tr.origname[3:].split('.')[0] for tr in superlocus.transcripts if tr.orig_annot]
		
		# first simple case - no original annotations, just give everyone novel names and call it a day
		if len(origtranscripts) == 0:
			for win in winners:
				if win.coding == True:
					win.addName("novel_coding_"+str(coding_stats[3]))
					coding_stats[3] += 1
				else:
					win.addName("novel_noncoding_"+str(noncoding_stats[3]))
					noncoding_stats[3] += 1				
		else:
			# second simple case - winner is an original annotation, keep its name (stripped of leading an_ and transcript #)
			remwinners = []; assigned = []
	#		print("First check - is this an original annot?")
			for win in winners:
				if win.orig_annot:
					newname = win.origname[3:].split('.')[0]
	#				print("Is original annot, assign name:",newname)
	#				print("Already assigned:",assigned)
					if newname in assigned:
						print("Warning: name",newname,"already assigned (this is either an error in the original annotations or a bug in the program)")
						print("Tried to assign name",newname,"to",win.origname)
						print("Current contents of assigned:")
						print(assigned)
					# note the addName function handles this
					win.addName(newname)
					assigned.append(newname)
					if win.coding == True:
						coding_stats[0] += 1
					else:
						noncoding_stats[0] += 1
				else:
					remwinners.append(win)
			
			winners = remwinners
			
			# for remainder, compare all winners to original annots to pick best name
			if len(winners) > 0:	
	#			print("Winners:",[a.origname for a in winners])
	#			print("Original:",[a.origname for a in origtranscripts])
		
				compscores = []; remwinners = []
				for win in winners:
					compscoretmp = []
					for origtr in origtranscripts:	
						comp = win.compareTranscripts(origtr)
	#					print("Comparison of",win.origname,"to",origtr.origname,"is:",comp)
						compscoretmp.append(comp[0])
								
					# deal with easy cases here - should be most of them
					if (win.coding == True and min(compscoretmp) >= 7) or (win.coding == False and min(compscoretmp) >= 4):
						# little to no overlap with any original annotation -> novel
	#					print("No strong overlap with any existing annotations, assign novel name")
						if win.coding == True:
							win.addName("novel_coding_"+str(coding_stats[3]))
							coding_stats[3] += 1
						else:
							win.addName("novel_noncoding_"+str(noncoding_stats[3]))
							noncoding_stats[3] += 1				
					elif (win.coding == True and min(compscoretmp) <= 3) or (win.coding == False and min(compscoretmp) <= 1):
						# perfect or nearly perfect match -> match
						matchidx = compscoretmp.index(min(compscoretmp))
						newname = origtranscripts[matchidx].origname[3:].split('.')[0]
	#					print("Perfect or nearly perfect match to existing annotation",newname,", assigning name ->",newname)
						if newname in assigned:
							print("Warning: name",newname,"already assigned (this is either an error in the original annotations or a bug in the program)")
							print("Tried to assign name",newname,"to",win.origname)
							print("Current contents of assigned:")
							print(assigned)
						win.addName(newname)
						assigned.append(newname)
						if win.coding == True:
							coding_stats[coding_mapping[min(compscoretmp)]] += 1
						else:
							noncoding_stats[coding_mapping[min(compscoretmp)]] += 1
					else:
	#					print("Still needs a name:",win.origname)
						remwinners.append(win)
						compscores.append(compscoretmp)
					
					winners = remwinners

			# remaining cases
			if len(winners) > 0:				
				# store comparison scores of all combinations of winners and original transcripts
				# bestnames[winner] = list of origtranscripts that overlap winner
				bestnames = []; matchidxs = []
				for i in range(len(winners)):
					bestname = []; matchidx = []
					for j in range(len(origtranscripts)):
						if (winners[i].coding == True and compscores[i][j] < 8) or (winners[i].coding == False and compscores[i][j] < 5):
							stripname = origtranscripts[j].origname[3:].split('.')[0]
							if stripname not in bestname:
								bestname.append(stripname)
								matchidx.append(j)
					
					bestnames.append(bestname)
					matchidxs.append(matchidx)

				# if two (or more) original annotations were merged, then multiple original annotations will overlap the winner
				if len(winners) > 1 or len(origtranscripts) > 1:
	#				print("Multiple winner and/or orig transcripts left, check for fusions/splits")						
					for i in range(len(winners)):
						if len(bestnames[i]) > 1:
	#						print(winners[i].origname,"is a possible merge of these original annots:",bestnames[i])					
							
							# find which of the original transcripts are fully or nearly fully contained in the novel one
							containedidx = []
							for j in range(len(bestnames[i])):
	#							print("Winner",winners[i].origname,"has exons:",winners[i].exons)
	#							print("Comparing to:",origtranscripts[matchidxs[i][j]].origname,"with exons",origtranscripts[matchidxs[i][j]].exons)
							
								if winners[i].overlap_by_exon(origtranscripts[matchidxs[i][j]])/origtranscripts[matchidxs[i][j]].totlen > 0.9:
	#								print(origtranscripts[matchidxs[i][j]].origname,"is fully contained")
									containedidx.append(j)
		
							#  make sure none of the contained transcripts fully overlap each other (shouldn't happen... but original annots may have errors)
							containednames = []; contain = []
							for z in range(len(containedidx)):
	#							print("Checking if contained gene",bestnames[i][containedidx[z]],"is fully contained by others...")
								orig_contained = False
								for h in range(len(containedidx)):
									if z != h:
										if origtranscripts[matchidxs[i][containedidx[z]]].overlap_by_exon(origtranscripts[matchidxs[i][containedidx[h]]]) / origtranscripts[matchidxs[i][containedidx[z]]].totlen > 0.9:
	#										print("Original transcript",bestnames[i][containedidx[z]],"is fully contained in",bestnames[i][containedidx[h]])
											orig_contained = True
								if orig_contained == False:
									contain += origtranscripts[matchidxs[i][containedidx[z]]].exons
									containednames.append(bestnames[i][containedidx[z]])
								
									
	#						print(winners[i].origname,"fully or almost fully contains",containednames)
							
							# test whether these original transcripts combined together account for more than 80% of the novel transcript;
							# if yes it's a valid merge, if no then continue
							if len(containednames) > 1:
								fracexplained = winners[i].overlap_by_exon2(contain) / winners[i].totlen							
								if fracexplained > 0.8:
									merged += 1
									winners[i].type = "Merged_annot"
								
									# if only two genes, give gene name that is the two combined
									if len(containednames) == 2:
										newname = "merge"
										for nn in containednames:
											newname += "_" + nn
											if nn in assigned:
												print("Warning: name",nn,"from merge already in assigned")
											assigned.append(nn)
										winners[i].addName(newname)			
	#									print("(Merge) Assigning name",winners[i].finalname,"to winner",winners[i].origname)
									elif len(containednames) > 2:
										# else give it 'novel gene' name
										winners[i].addName("novel_merged_"+str(merged))							
	#									print("(Merge) Assigning name",winners[i].finalname,"to winner",winners[i].origname)
								
					# if an original annotation was split into two (or more) winners, then the same original annotation will overlap multiple winners
					# can only happen if two or more winners remain
					if len(winners) >= 2:
	#					print("Looking into split of these:",[w.origname for w in winners])
	#					print(bestnames)
	#					print("Original annots at this locus:")
	#					print(origtranscriptnames)
						
						for orignames in set([j for sub in bestnames for j in sub]):
	#						print("Could original annot",orignames,"be split?")
	#						print([nn for nn in bestnames if orignames in nn])
						
							if len([nn for nn in bestnames if orignames in nn]) >= 2:
								# this original annot overlaps multiple winners, identify them
								idx = origtranscriptnames.index(orignames)
	#							print("Original annot",orignames,"overlaps multiple winners (index",idx,")")
								
								# get the winners that are fully or nearly fully contained in the original annot and add together their exons
								# to see how much they account for the full original annot
								contain = []; overlapidx = []
								for i in range(len(bestnames)):
									if winners[i].overlap_by_exon(origtranscripts[idx])/winners[i].totlen > 0.9:
	#									print("Winner",winners[i].origname,"is fully contained within",origtranscriptnames[idx])
										contain += winners[i].exons
										overlapidx.append(i)
								
	#							print("Index of strongly overlapping winners:",overlapidx)
	#							print("Full list of winners:",[x.origname for x in winners])
								
								# test whether these 2+ winners combined together account for more than 80% of the original transcript;
								# if yes it's a valid split, if no then continue
								fracexplained = origtranscripts[idx].overlap_by_exon2(contain) / origtranscripts[idx].totlen									
								if fracexplained > 0.8 and len(overlapidx) > 1:
									ssplit += 1
									nn = 1
									if orignames in assigned:
										print("Warning: name",orignames,"from split already in assigned")
									assigned.append(orignames)
									for oidx in overlapidx:
										winners[oidx].type = "Split_annot"
										newname = orignames + "_s" + str(nn)
										winners[oidx].addName(newname)
										print("(Split) Assigning name",newname,"to winner",winners[oidx].origname)
										nn += 1

				# Get what's left
				todo = []
				for i in range(len(winners)):
					if winners[i].finalname is None:
						todo.append(i)
						
				# Everything left is not a split or a merge, so just assign best if not yet assigned
				# and they're not too different, or assign novel name
				if len(todo) > 0:
	#				print("Multiple winners left:")
	#				print("All winners:")
	#				print([a.origname for a in winners])
	#				print("Still need final name:")
	#				print([winners[i].origname for i in todo])
	#				print("Compare to original transcripts:")
	#				print([a.origname for a in origtranscripts])
	#				print(compscores)
	#				print("Filtered compscores:")
	#				print([compscores[i] for i in todo])
	#				print("Already assigned:",assigned)
				
					for i in range(len(todo)):
	#					print("Finding best name for",winners[todo[i]].origname)
	#					print("Best score is",min(compscores[todo[i]]))
	#					print("Orig matching best score is",origtranscripts[compscores[todo[i]].index(min(compscores[todo[i]]))].origname[3:].split('.')[0])
					
						bestname = origtranscripts[compscores[todo[i]].index(min(compscores[todo[i]]))].origname[3:].split('.')[0]	
						
						if winners[todo[i]].coding == True and min(compscores[todo[i]]) <= 6 and bestname not in assigned:
	#						print("Assigning name",bestname)
							winners[todo[i]].addName(bestname)
							coding_stats[coding_mapping[min(compscores[todo[i]])]] += 1
						elif winners[todo[i]].coding == False and min(compscores[todo[i]]) <= 3 and bestname not in assigned:
	#						print("Assigning name",bestname)
							winners[todo[i]].addName(bestname)
							noncoding_stats[noncoding_mapping[min(compscores[todo[i]])]] += 1
						else:
							if winners[todo[i]].coding == True:
	#							print("Assigning name novel_coding_"+str(coding_stats[3]))
								winners[todo[i]].addName("novel_coding_"+str(coding_stats[3]))
								coding_stats[3] += 1
							else:
	#							print("Assigning name novel_noncoding_"+str(coding_stats[3]))
								winners[todo[i]].addName("novel_noncoding_"+str(noncoding_stats[3]))
								noncoding_stats[3] += 1				
								

# filters transcript list to try to improve recursion levels, by dropping transcripts that
# are very similar to others in the group but score more poorly
def filterTranscripts(ss, target = 30):

	todrop = len(ss.transcripts) - target
#	print("Trying to filter set of",len(ss.transcripts),"transcripts:")
#	print([tt.origname for tt in ss.transcripts])

	# first, assign a score to all transcripts in superlocus ss
	idxcombos = [str(i) for i in range(len(ss.transcripts))]
	scorearr,maxidx,maxscore = getScores(ss, idxcombos, outputstats = False)
	
#	print("Transcripts:",[x.origname for x in ss.transcripts])
#	print("Scores:",scorearr)
	
	# now get all sets of strongly overlapping transcripts
	visited = [0]*len(ss.transcripts); overlapsets = []
	for i in range(len(ss.transcripts)):
		if visited[i] == 0 and i < len(ss.transcripts):
	#		print("Finding all transcripts that strongly overlap with",ss.transcripts[i].origname,", idx =",i)
			visited[i] = 1
			# find all unvisited transcripts similar to this one
			overlset = [i]
			for j in range(i+1,len(ss.transcripts)):
				if visited[j] == 0:
					overl = ss.transcripts[i].overlap_by_exon(ss.transcripts[j])
					overlfrac1 = min(overl / ss.transcripts[i].totlen, overl / ss.transcripts[j].totlen)
					overlfrac2 = max(overl / ss.transcripts[i].totlen, overl / ss.transcripts[j].totlen)
	#				print("Overlap between",ss.transcripts[i].origname,"and",ss.transcripts[j].origname,",",overlfrac1,overlfrac2)
					if overlfrac1 > 0.8 or (overlfrac2 > 0.99 and overlfrac1 > 0.6):
	#					print(" -> ",ss.transcripts[j].origname,"overlaps (idx =",j,")")
						visited[j] = 1
						overlset.append(j)
	#		print("Final overlap set:",overlset)
			overlapsets.append(overlset)
		
	# now for each set, starting from the biggest one, drop bottom x transcripts until
	# either only one remains (move on to next set) or todrop reached
	overlapsets.sort(key=len)
	overlapsets = overlapsets[::-1]
#	print("Overlapsets:",overlapsets)
	
	dropped = 0; todropidx = []
	for overlset in overlapsets:
#		print("Overlap set:",overlset)
		scores = [scorearr[i] for i in overlset]
#		print("Scores:",scores)
		sortidx = [sidx[0] for sidx in sorted(enumerate(scores),key=lambda i:i[1])]
		sortscores = [scores[i] for i in sortidx]
		sortset = [overlset[i] for i in sortidx]
		for i in range(len(sortset)-1):
			if dropped >= todrop:
				break
			else:
				dropped += 1
				todropidx.append(sortset[i])
		if dropped >= todrop:
			break

	# if we weren't able to drop enough this way, instead start deleting the lowest scoring overall instead of by overlapping set
	if len(todropidx) < todrop:
	#	print("Dropping lowest scoring overall:")
	#	print("Already dropping:",todropidx)
	#	print(scorearr)
	#	print([tt.origname for tt in ss.transcripts])
		sortidx = [sidx[0] for sidx in sorted(enumerate(scorearr),key=lambda i:i[1])]
		sortscores = [scorearr[i] for i in sortidx]
		for i in range(len(sortidx)-1):
			if dropped >= todrop:
				break
			else:
				dropped += 1
				todropidx.append(sortidx[i])
		
	# drop these transcripts
#	print("To drop:",todropidx)		
	for idx in sorted(todropidx, reverse = True):								# drop in reverse order so indexes aren't thrown off
#		print("Dropping index",idx,", name",ss.transcripts[idx].origname)
		del ss.transcripts[idx]
	
	# check whether deleting this many transcripts improved getValidIntervalsFuzzy() run
	tt = [x.mRNA for x in ss.transcripts]
	idxcombos = []
	getValidIntervalsFuzzy(tt, 0, -1, "", idxcombos, True)
			
	if len(idxcombos) == args.max_recursion + 1:
		if target - 5 > 0:
			print("Initial filtering not sufficient, still hitting recursion threshold; try reducing target to",target - 3)
			res = filterTranscripts(ss, target-3)
			return(res)
		else:
			print("Failed to hit recursion target and can't drop more transcripts (this really shouldn't happen); only",len(ss.transcripts),"remain")
			return(ss,False)
	else:
		print("Recursion successful with",target,"transcripts remaining")
		return(ss,True)


def mergeIntervals(intervals):
	intervals.sort(key = lambda ex: ex[0])
	mergedintervals = [intervals[0]]
	for cur in intervals:
		prev = mergedintervals[-1]
		if cur[0] <= prev[1]:
			prev[1] = max(prev[1], cur[1])
		else:
			mergedintervals.append(cur)	
	return mergedintervals		
	

def getOutside(exonlist,boundaries,chr,strand):
	# given a list of exons contained within an outer boundary, returns the average expression
	# over regions in the boundary not covered by the exons
#	print("Running getOutside() on",exonlist,"w/ boundaries",boundaries,"on chr",chr)
	
	exonlist = [list(ex) for ex in exonlist]			# convert to list of lists if not already and sort
	exonlist.sort(key = lambda ex: ex[0])
	mergedexonlist = mergeIntervals(exonlist)		
	
	# subtract exonlist from superlocus boundaries to get all intervals not covered by exons
	outsideidx = []
	left = boundaries[0]
	for ex in mergedexonlist:
		if left < ex[0]-1:
			outsideidx.append((left,ex[0]-1))
		left = ex[1]+1
	if left < boundaries[1]:
		outsideidx.append((left,boundaries[1]))
		
	# drop 'exons' that are too short:
	outsideidx = [ex for ex in outsideidx if ex[1]-ex[0]+1 >= 20]

	# get weighted average coverage over all outside intervals (weighted by exon length)
	if len(outsideidx) > 0:
		if strand == "+":
			otmp = [plustrack.stats(chr, ex[0], ex[1])[0] for ex in outsideidx]
		else:		
			otmp = [minustrack.stats(chr, ex[0], ex[1])[0] for ex in outsideidx]

		# get weighted avg coverage (weight by exon length)
		tmp1 = 0; tmp2 = 0
		for i in range(len(otmp)):
			tmp1 += otmp[i]*(outsideidx[i][1]-outsideidx[i][0]+1)
			tmp2 += (outsideidx[i][1]-outsideidx[i][0]+1)
		return tmp1,tmp2		# return sum and length separately for weighted averaging purposes
	else:
		return 0,0


def getScores(superlocus, idxcombos, outputstats = False):
	
	# for each combination in idxcombos, get an average overall score combining (across all transcripts in the combo:
	# CDS length (longer is better), coding >> noncoding if both are present
	# Avg expression (more is better)
	# Number of junction reads supported
	# Picked by Mikado
	# Original annot is also given some weight	
	# Outside expression (not captured by annots)
	cdslen = []; totlen = []; totcov = []; avgcov = []; stdevcov = []; juncreads = []; utrdist = []; isorig = []; ispicked = []; 

	for idx in idxcombos:
#		print("Getting score for",idx)
		idxarr = [int(x) for x in idx.split(';')]
		c = 0; t = 0; tc = []; a = []; s = []; j = 0; uu = 0; og = 0; pi = 0
		for ii in idxarr:
			# if a CDS is crappy (borderline passes cutoffs) then pretend it's not there for scoring (don't want to accidentally give huge priority to a tiny crappy protein vs. a nice ncRNA)
			if not superlocus.transcripts[ii].orig_annot and not superlocus.transcripts[ii].forcedCDS and (superlocus.transcripts[ii].totCDS < 1.5*args.min_prot_len or superlocus.transcripts[ii].totCDS/superlocus.transcripts[ii].totlen < args.min_orf_frac*1.5):
				c += 0
			else:
				c += superlocus.transcripts[ii].totCDS
			t += superlocus.transcripts[ii].totlen
			tc.append(superlocus.transcripts[ii].tot_cov)
			a.append(superlocus.transcripts[ii].avg_cov)
			s.append(superlocus.transcripts[ii].var_all)
			j += superlocus.transcripts[ii].tot_supp_juncs
			if superlocus.transcripts[ii].coding:
				uu += superlocus.transcripts[ii].UTR_dist
			if superlocus.transcripts[ii].picked == True:
				pi = 1
			if superlocus.transcripts[ii].orig_annot == True:
				og = 1
		
		cdslen.append(c)
		totlen.append(t)
		avgcov.append(sum(a)/len(a))
		totcov.append(sum(tc))
		stdevcov.append(sum(s)/len(s))
		juncreads.append(j)
		utrdist.append(uu)
		isorig.append(og)
		ispicked.append(pi)
		
	# add one more statistic: read coverage remaining inside superlocus but outside of transcript(s)
	# in each combo, to try to prioritize combos that capture the majority of reads in a locus
	# (sums across all regions -not- occupied by exons of transcript(s) in combo)
	outside = []; outside_int = []
	for idx in idxcombos:
		idxarr = [int(x) for x in idx.split(';')]
	
		# get coverage over introns ('internal' outside coverage)
		tot_outside_internal = 0; len_outside_internal = 0
		for ii in idxarr:
			res = getOutside(superlocus.transcripts[ii].exons,superlocus.transcripts[ii].mRNA,superlocus.chr,superlocus.strand)
			tot_outside_internal += res[0]; len_outside_internal += res[1]
		if len_outside_internal == 0:
			outside_internal = 0
		else:
			outside_internal = tot_outside_internal/len_outside_internal
		
		# get coverage outside of gene bounds
		intervals = []
		for ii in idxarr:
			intervals.append(list(superlocus.transcripts[ii].mRNA))
		res = getOutside(intervals,superlocus.boundaries,superlocus.chr,superlocus.strand)
		if res[1] == 0:
			outside_external = 0
		else:
			outside_external = res[0] / res[1]
							
#		if "cl_Chr1.5556.0" in [superlocus.transcripts[x].origname for x in idxarr]:
#			print("Getting outside_internal and outside_external values for indexes",idxarr)
#			print("Gene(s):",[superlocus.transcripts[x].origname for x in idxarr])
#			print("External coverage:",outside_external)
#			print("Internal coverage:",outside_internal)
		
		outside.append(outside_external)
		outside_int.append(outside_internal)

	# nice table output summary for all transcript combos in a superlocus, showing all attributes used to evaluate
	if outputstats == True:
		print('index\tcdslen\ttotlen\ttotcov\tavgcov\tstdevcov\tjuncreads\tUTRdist\tisorig\tispicked\toutside_int\toutside_ext\ttranscript_names')
		ii = 0
		for idx in idxcombos:
			idxarr = [int(x) for x in idx.split(';')]
			nn = ','.join([superlocus.transcripts[j].origname for j in idxarr])
			print(idx,'\t',cdslen[ii],'\t',totlen[ii],'\t',totcov[ii],'\t',avgcov[ii],'\t',stdevcov[ii],'\t',juncreads[ii],'\t',utrdist[ii],'\t',isorig[ii],'\t',ispicked[ii],'\t',outside_int[ii],'\t',outside[ii],'\t',nn)
			ii+=1
			if ii > 500:
				print("Reached printing limit of 500, additional rows not shown")
				break

	# scale everything to 0-100 and round to nearest whole for readibility		
	# for these attributes, higher is better, so keep as-is
	# also, if a CDS is particularly short or a small fraction of the gene (but still valid),
	# treat as noncoding here as it's not a very 'good' CDS
	if max(cdslen) == 0:
		cdslenscaled = cdslen		# if all zeros, scaling sets them all to 100 -> bad
	else:
		cdslenscaled = np.array(cdslen)
		cdslenscaled = np.round(np.interp(cdslenscaled, (0,cdslenscaled.max()), (0,100)),0).tolist()
	
	totcovscaled = np.array(totcov)
	if totcovscaled.max() < args.min_denovo_expr:
		totcovscaled = [0]*len(totcovscaled)
	else:		
		totcovscaled = np.round(np.interp(totcovscaled, (0,totcovscaled.max()), (0,100)),0).tolist()

	# note - for this and outside, if the highest value is too low (<= min_denovo_expr) then just set to 0 all around
	avgcovscaled = np.array(avgcov)
	if avgcovscaled.max() < args.min_denovo_expr:
		avgcovscaled = [0]*len(avgcovscaled)
	else:		
		avgcovscaled = np.round(np.interp(avgcovscaled, (0,avgcovscaled.max()), (0,100)),0).tolist()

	# if nearly no junction reads detected, just ignore
	juncreadsscaled = np.array(juncreads)
	if juncreadsscaled.max() < 3:
		juncreadsscaled = [0]*len(juncreadsscaled)
	else:
		juncreadsscaled = np.round(np.interp(juncreadsscaled, (0,juncreadsscaled.max()), (0,100)),0).tolist()
	
	# for these attributes, lower is better, so after scaling, flip it around
	totlenscaled = np.array(totlen)
	totlenscaled = (100-np.round(np.interp(totlenscaled, (0,totlenscaled.max()), (0,100)),0)).tolist()
		
	stdevcovscaled = np.array(stdevcov)
	if stdevcovscaled.max() < args.min_denovo_expr:
		stdevcovscaled = [0]*len(stdevcovscaled)
	else:		
		stdevcovscaled = (100-np.round(np.interp(stdevcovscaled, (0,stdevcovscaled.max()), (0,100)),0)).tolist()

	utrdistscaled = np.array(utrdist)
	utrdistscaled = (100-np.round(np.interp(utrdistscaled, (0,utrdistscaled.max()), (0,100)),0)).tolist()

	# note - for this and avgcov, if the highest value is 0 then just set to 0 all around
	outsidescaled = np.array(outside)
	if outsidescaled.max() < args.min_denovo_expr:
		outsidescaled = [0]*len(outsidescaled)
	else:		
		outsidescaled = (100-np.round(np.interp(outsidescaled, (0,outsidescaled.max()), (0,100)),0)).tolist()

	# note - for this and avgcov, if the highest value is 0 then just set to 0 all around
	outsideintscaled = np.array(outside_int)
	if outsideintscaled.max() < args.min_denovo_expr:
		outsideintscaled = [0]*len(outsideintscaled)
	else:		
		outsideintscaled = (100-np.round(np.interp(outsideintscaled, (0,outsideintscaled.max()), (0,100)),0)).tolist()

	# for both isorig and ispicked, set to 100 if 1, 0 otherwise
	isorig = [x*100 for x in isorig]
	ispicked = [x*100 for x in ispicked]		
							
	# calculate overall score according to user-provided weights
	maxscore = 0; maxidx = -1; ii = 0; scorearr = []
	for idx in idxcombos:
		score = cdslenscaled[ii]*args.weight_cdslen + totlenscaled[ii]*args.weight_totlen + totcovscaled[ii]*args.weight_totcov + avgcovscaled[ii]*args.weight_avgcov + stdevcovscaled[ii]*args.weight_stdevcov + juncreadsscaled[ii]*args.weight_juncreads + utrdistscaled[ii]*args.weight_utrdist + isorig[ii]*args.weight_isorig + ispicked[ii]*args.weight_ispicked + outsideintscaled[ii]*args.weight_outside_internal + outsidescaled[ii]*args.weight_outside
		scorearr.append(score)
		if score > maxscore:
			maxscore = score
			maxidx = ii
		elif score == maxscore:
			# break ties by favoring existing annotation, then by CDS len, then by total len, then random
			if isorig[ii] == 1 and isorig[maxidx] == 0:
				maxidx = ii
			elif cdslen[ii] > cdslen[maxidx]:
				maxidx = ii
			elif totlen[ii] > totlen[maxidx]:
				maxidx = ii
		ii += 1

	if outputstats == True:
		print('')
		print('index\tcdslen\ttotlen\ttotcov\tavgcov\tstdevcov\tjuncreads\tutrdist\tisorig\tispicked\toutside_int\toutside_ext\tscore\ttranscript_names')
		ii = 0
		for idx in idxcombos:
			idxarr = [int(x) for x in idx.split(';')]
			nn = ','.join([superlocus.transcripts[j].origname for j in idxarr])
			print(idx,'\t',cdslenscaled[ii],'\t',totlenscaled[ii],'\t',totcovscaled[ii],'\t',avgcovscaled[ii],'\t',stdevcovscaled[ii],'\t',juncreadsscaled[ii],'\t',utrdistscaled[ii],'\t',isorig[ii],'\t',ispicked[ii],'\t',outsideintscaled[ii],'\t',outsidescaled[ii],'\t',scorearr[ii],'\t',nn)
			ii+=1
			if ii > 500:
				print("Reached printing limit of 500, additional rows not shown")
				break

		print("Max score was:",maxscore,"for combo",idxcombos[maxidx])
		
	# found best combo! return full list of scores as well as the index and value of the max (breaking ties as shown above)
	return scorearr,maxidx,maxscore


# identifies, renames and outputs best transcripts for each superlocus
def pickBestTranscript(superlocus):

	global usednames
	global coding_stats
	global noncoding_stats

	if args.printstats != "" and args.printstats in [a.origname for a in superlocus.transcripts]:
		print("Getting best transcript(s) for superlocus:")
		superlocus.printInfo()
		for tt in superlocus.transcripts:
			tt.printTranscript()
		
	# if there's just one annotation, print it unless it is NOT an original annotation
	# and has low average coverage
	if len(superlocus.transcripts) == 1:
		superlocus.transcripts[0].winner = True
		if args.printstats != "" and args.printstats in [a.origname for a in superlocus.transcripts]:
			print("This superlocus has only one locus, so it automatically wins!")
	else:
		# else, find best combination of transcript for this loci
	
		# (1) First, get all possible combinations of nonoverlapping transcripts
		# TODO: sub in intersection matrix, to allow noncoding things in introns (rare in plants tho)
	
		# make intersection matrix - True/False whether or not two transcripts are overlapping
		# TODO for modified getValidIntervals function, not currently used
#COM		intersect_matrix = [[False for x in range(len(tt))] for y in range(len(tt))]
#COM		for i in range(len(superlocus.transcripts)):
#COM			for j in range(i+1,len(superlocus.transcripts)):
#COM#				print("Overlap between:")
#COM#				superlocus.transcripts[i].printTranscript()
#COM#				superlocus.transcripts[j].printTranscript()
#COM				
#COM				overlap = superlocus.transcripts[i].overlap_by_exon(superlocus.transcripts[j])
#COM				if overlap <= 40 and max(overlap/superlocus.transcripts[i].totlen, overlap/superlocus.transcripts[j].totlen) < 0.2:
#COM					intersect_matrix[i][j] = True
#COM					intersect_matrix[j][i] = True
		
		tt = [x.mRNA for x in superlocus.transcripts]		# get bounds for all transcripts
				
		idxcombos = []
		getValidIntervalsFuzzy(tt, 0, -1, "", idxcombos, True)
				
		if len(idxcombos) == args.max_recursion + 1:
			print("Max level of recursion",args.max_recursion,"hit in getValidIntervalsFuzzy() for superlocus",superlocus.chr+":"+str(superlocus.boundaries[0])+"-"+str(superlocus.boundaries[1]))
			print("Reducing complexity at this locus by dropping poor quality transcripts...")
	#		print("Original superlocus has",len(superlocus.transcripts),"transcripts")
			if len(superlocus.transcripts) > 3:
				superlocus,result = filterTranscripts(superlocus, len(superlocus.transcripts)-3)
			else:
				print("This locus already has few transcripts, can't drop more (and the recursion shouldn't be so high with so few transcripts...)")
				result = False
	#		print("Post filtering, superlocus has",len(superlocus.transcripts),"transcripts")
			
			if result == True:
	#			print("Recursion now successful with",len(superlocus.transcripts),"transcripts")			
				# repeat recursion to get intervals
				tt = [x.mRNA for x in superlocus.transcripts]
				idxcombos = []
				getValidIntervalsFuzzy(tt, 0, -1, "", idxcombos, True)				
			else:
				# recursion failed; don't consider any gene combos, just output a list of single indexes
				print("Warning: could not get transcript set to run under recursion limit, calls at locus",superlocus.chr+":"+str(superlocus.boundaries[0])+"-"+str(superlocus.boundaries[1]),"may be suboptimal")
				idxcombos = []
				for i in range(len(superlocus.transcripts)):
					idxcombos.append(str(i))
				
			
		# for each combination in idxcombos, get overall score using getScores()
		outputstats = False
		if args.printstats != "" and args.printstats in [a.origname for a in superlocus.transcripts]:
			outputstats = True
		
		# if there are 4 or more idxcombos, then divide and conquer (this is to eliminate an effect
		# where very bad combos set the scale of the normalized attributes used for scoring, and make
		# better combos all look the same so that that attribute doesn't matter in selection of final winner
		scorearr,maxidx,maxscore = getScores(superlocus, idxcombos, outputstats)
		
		while len(idxcombos) >= 4:
			# drop the bottom half worst scoring idxcombos, then repeat
			if args.printstats != "" and args.printstats in [a.origname for a in superlocus.transcripts]:
				print("Refining calls by dropping bottom",math.floor(len(scorearr)/2),"out of",len(scorearr),"index combos and recalculating scores...")
			
			sortidx = [sidx[0] for sidx in sorted(enumerate(scorearr),key=lambda i:i[1])]
			tokeep = sortidx[math.floor(len(sortidx)/2):len(sortidx)]
			idxcombos = [idxcombos[i] for i in tokeep]			
			scorearr,maxidx,maxscore = getScores(superlocus, idxcombos, outputstats)
										
		
		# found best combo! now update names for each and print to file
		idxarr = [int(x) for x in idxcombos[maxidx].split(';')]
		for ii in idxarr:
			superlocus.transcripts[ii].winner = True
		
		# additionally, if there are any original loci in this superlocus that were not picked and do
		# NOT overlap any of the picked annotations (generally just means they have no reads in the input data)
		# also select those here
		for tt in superlocus.transcripts:
			if tt.winner == False and tt.orig_annot == True:
	#			print(tt.origname,"is an orig transcript but was not a winner, check if it overlaps any of the winners")
				maxoverlap = 0
				for win in [w for w in superlocus.transcripts if w.winner == True]:
	#				print("Overlap between",tt.origname,"and",win.origname,"is",tt.overlap_by_exon(win))
					overlap = tt.overlap_by_exon(win)
					if overlap > maxoverlap:
						maxoverlap = overlap
	#			print("Max overlap between",tt.origname,"and any winner is",maxoverlap)
							
				# if transcripts don't overlap and tt is original but not winner, make a winner
				if maxoverlap <= args.fuzziness:	
					tt.winner = True
				
	#	if args.printstats != "" and args.printstats in [a.origname for a in superlocus.transcripts]:		
	#		sys.exit(1)
	
	# Winners have been decided, now output to file
#	print("Winners for superlocus at",superlocus.chr,":",superlocus.boundaries[0],"-",superlocus.boundaries[1],"have been assigned, outputting to file")						
	outputSuperlocus(superlocus)
		
		
		
		
		

#-------------------------------------------------------------
# MAIN CODE
#-------------------------------------------------------------

# Check all input files exist
toopen = [args.pick, args.raw, args.junc, args.genome, args.plustrack, args.minustrack]
toopenname = ["pick", "raw", "junc", "genome", "plustrack", "minustrack"]
for ff,ss in zip(toopen,toopenname):
	if ff is None:
		print("Error: option --"+ss,"is required")
		sys.exit(1)
	elif not os.path.exists(ff):
		print("Error: can't open file",ff)
		sys.exit(1)

# Print options
print("")
print("Running mikado_refine v1.0		by Colette L. Picard, 01/11/2021")
print("-------------------------")
print("Required input files:")
print("Picked mikado transcripts:",args.pick)
print("All mikado transcripts:",args.raw)
print("Portcullis junctions file:",args.junc)
print("Genome fasta file:",args.genome)
print("BigWig coverage track (+ strand):",args.plustrack)
print("BigWig coverage track (- strand):",args.minustrack)
print("---------")
print("Transcript selection options:")
print("Criteria for valid transcripts:")
print(" - Maximum total transcript length (including introns):", args.maxlen)
print(" - Minimum fraction exonic (spliced mRNA length / unspliced mRNA length):", args.min_exon_ratio)
print("ORF detection criteria:")
print(" - Minimum total CDS length:",args.min_prot_len)
print(" - Minimum percent CDS (CDS length / spliced mRNA length):",args.min_orf_frac)
print("---------")
print("Other options:")
print("Prefix for output files:",args.outprefix)
print("-------------------------")

print("")
print("Preprocessing steps:")

# Open genome and bigWig files
print("Loading genome and bigWig files...")
genomeseq = pysam.FastaFile(args.genome)
plustrack = pyBigWig.open(args.plustrack)
minustrack = pyBigWig.open(args.minustrack)

# Open and read in the portcullis junctions file
print("Loading Portcullis junction file...")
try:
	f = open(args.junc, 'r')
except IOError as e:
	print(e)
	print('Could not open input file',args.junc)
	sys.exit(2)

juncs = {}			# store junctions in 2-lvl dict by chr+strand then by start -> store (end, depth) tuple

line = f.readline(); lnum = 0
line = f.readline()					# ignore first line (BED file header)
while line:
	lnum += 1
#	print(line.strip())
	
	ll = line.strip().split('\t')
	chr = ll[0]
	strand = ll[5]
	depth = int(float(ll[4]))
	jstart = int(ll[6])
	jend = int(ll[7]) + 1 			# with this indexing, LHS of junction == RHS of exon 1, RHS junc == LHS of exon 2
	
	if chr+strand not in juncs:
		juncs[chr+strand] = {}
	if jstart in juncs[chr+strand]:
		juncs[chr+strand][jstart].append((jend,depth))
	else:
		juncs[chr+strand][jstart] = [(jend,depth)]
	
	line = f.readline()
f.close()


# Open and read in the picked mikado GFF3 file, store in dict for fast lookup keyed
# to overall transcript boundary
print("Reading in mikado pick output...")
try:
	f = open(args.pick, 'r')
except IOError as e:
	print(e)
	print('Could not open input file',args.pick)
	sys.exit(2)

mikado_picked = {}			# store all mikado pick data in dict indexed by chr+strand, then by transcript bounds (transcript.mRNA)

tmp_transcript = None
cur_transcript = ""
line = f.readline()
while line:
	if line.strip() and line[0] != '#':
		ll = line.strip().split('\t')	
		if ll[2] == "exon":
			transcript_id = ll[8].split('Parent=')[-1]
			istart = int(ll[3])
			iend = int(ll[4])		
		
			if transcript_id != cur_transcript:
				if tmp_transcript is not None:
					# current transcript is complete, finish it by updating mRNA bounds, then add to dict
					tmp_transcript.finishPartialEntry()
					if tmp_transcript.chr+tmp_transcript.strand not in mikado_picked:
						mikado_picked[tmp_transcript.chr+tmp_transcript.strand] = {}
					if tmp_transcript.mRNA in mikado_picked[tmp_transcript.chr+tmp_transcript.strand]:
						mikado_picked[tmp_transcript.chr+tmp_transcript.strand][tmp_transcript.mRNA].append(tmp_transcript)
					else:
						mikado_picked[tmp_transcript.chr+tmp_transcript.strand][tmp_transcript.mRNA] = [tmp_transcript]
			
				# create new transcript object from the new line
				tmp_transcript = Transcript(transcript_id, ll[0], ll[6], "picked")
				tmp_transcript.addLine("exon", istart, iend)
				cur_transcript = transcript_id		
			else:
				tmp_transcript.addLine("exon", istart, iend)	
	line = f.readline()

# done reading through file, add the last one
if tmp_transcript is not None:
	tmp_transcript.finishPartialEntry()
	if tmp_transcript.chr+tmp_transcript.strand not in mikado_picked:
		mikado_picked[tmp_transcript.chr+tmp_transcript.strand] = {}
	if tmp_transcript.mRNA in mikado_picked[tmp_transcript.chr+tmp_transcript.strand]:
		mikado_picked[tmp_transcript.chr+tmp_transcript.strand][tmp_transcript.mRNA].append(tmp_transcript)
	else:
		mikado_picked[tmp_transcript.chr+tmp_transcript.strand][tmp_transcript.mRNA] = [tmp_transcript]
f.close()

	
# Done with preprocessing, open the main mikado raw GTF file and process locus-by-locus
print("Done with preprocessing")
print("")
print("Going through each transcript in --raw file and picking best transcript(s) for each locus...")

try:
	f = open(args.raw, 'r')
except IOError as e:
	print(e)
	print('Could not open input file',args.raw)
	sys.exit(2)

try:
	outf = open(args.outprefix+'_final.gtf', 'w')
	outf_noAS = open(args.outprefix+'_final_noAS.gtf', 'w')
except IOError as e:
	print(e)
	print('Could not create output file',args.outprefix+'_final.gtf')
	sys.exit(2)

processed_genes = []
cur_transcript = ""
tmp_transcript = None
superlocusplus = None
superlocusminus = None
num_lines = sum(1 for line in open(args.raw))

# while loci are output as it goes, save superloci (with only winner transcripts remaining)
# so that at the end I can output a version without major noncoding antisense annotations (for stranded RNA-seq)
# store as: all_winners[chr+strand][leftidx / 100000] = list of transcripts in this region
all_winners = {}
chrlist = []

def save_winners(ss):
	global all_winners
	global chrlist

	winners = [tt for tt in ss.transcripts if tt.winner == True]
	if ss.chr+ss.strand not in all_winners:
		all_winners[ss.chr+ss.strand] = {}
	for tt in winners:
		leftl = math.floor(tt.mRNA[0]/100000)
		if leftl not in all_winners[ss.chr+ss.strand]:
			all_winners[ss.chr+ss.strand][leftl] = [tt]
		else:
			all_winners[ss.chr+ss.strand][leftl].append(tt)
			
	if ss.chr not in chrlist:
		chrlist.append(ss.chr)
		

line = f.readline(); lnum = 0
while line:
#	print(line.strip())
	lnum += 1

	if lnum % 1000 == 0:
		print("Processing line",lnum,"out of",num_lines,"("+str(round(lnum/num_lines*100,2))+"%)")
	
	if line.strip() and line[0] != '#':			# skip  header lines
		# parse 9th (last) field which has gene ID in it
		ll = line.strip().split('\t')
		transcript_id = ll[8].split('transcript_id ')[-1].split(';')[0][1:-1]
		chromosome = ll[0]
		source = ll[1]
		type = ll[2]
		start = int(ll[3])
		end = int(ll[4])
		strand = ll[6]

		# only process lines with these values in 3rd column
		if type in ["mRNA", "transcript", "exon", "CDS"]:
			
			# if this is a new transcript, finish current transcript and create new Transcript
			if transcript_id != cur_transcript:
			
				# proceed if this is a valid transcript
				if not tmp_transcript is None and tmp_transcript.isValid():
					tmp_transcript.finishFullEntry()
					
	#				if tmp_transcript.origname == "an_AT1G64390.1":
	#					tmp_transcript.printTranscript()
										
					# add to appropriate superlocus
					if tmp_transcript.strand == '+':
						if superlocusplus == None:			# first time only
							superlocusplus = Superlocus(tmp_transcript)
						elif tmp_transcript.overlapping2(superlocusplus) >= args.fuzziness:
							# add to current superlocus
	#						print("Adding",tmp_transcript.origname,"on strand",tmp_transcript.strand,"to + superlocus at",superlocusplus.boundaries,"on",superlocusplus.strand,"; overlap is",tmp_transcript.overlapping2(superlocusplus))
							superlocusplus.addTranscript(tmp_transcript)
						else:
							# current superlocus is complete, pick best transcripts from it, then create a new one with tmp_transcript
							superlocusplus.finishSuperlocus()
							pickBestTranscript(superlocusplus)
							save_winners(superlocusplus)										
							superlocusplus = Superlocus(tmp_transcript)
	#						print("Creating new + superlocus and adding",tmp_transcript.origname,"on strand",tmp_transcript.strand)
					else:
						if superlocusminus == None:			# first time only
							superlocusminus = Superlocus(tmp_transcript)
						elif tmp_transcript.overlapping2(superlocusminus) >= args.fuzziness:
							# add to current superlocus
	#						print("Adding",tmp_transcript.origname,"on strand",tmp_transcript.strand,"to - superlocus at",superlocusminus.boundaries,"on",superlocusminus.strand,"; overlap is",tmp_transcript.overlapping2(superlocusminus))
							superlocusminus.addTranscript(tmp_transcript)
						else:
							# current superlocus is complete, pick best transcripts from it, then create a new one with tmp_transcript
							superlocusminus.finishSuperlocus()
							pickBestTranscript(superlocusminus)
							save_winners(superlocusminus)
							superlocusminus = Superlocus(tmp_transcript)
	#						print("Creating new - superlocus and adding",tmp_transcript.origname,"on strand",tmp_transcript.strand)

				# create new transcript object for the new line
				tmp_transcript = Transcript(transcript_id, chromosome, strand, source)
				tmp_transcript.addLine(type, start, end)
				cur_transcript = transcript_id
			
			else:
				# not a new transcript, add this line to existing transcript
				tmp_transcript.addLine(type, start, end)
			
	line = f.readline()
	
f.close()

# done reading through, add the last entry and finish the last + and - superloci
if not tmp_transcript is None and tmp_transcript.isValid():
	if tmp_transcript.finished == False:
		tmp_transcript.finishFullEntry()
		
#	if tmp_transcript.origname == "an_AT1G68500.1":
#		tmp_transcript.printTranscript()
	
	if tmp_transcript.strand == '+':
		if superlocusplus == None:			# first time only
			superlocusplus = Superlocus(tmp_transcript)
		elif tmp_transcript.overlapping2(superlocusplus) >= args.fuzziness:
			superlocusplus.addTranscript(tmp_transcript)
		else:
			if superlocusplus is not None:
				superlocusplus.finishSuperlocus()
				pickBestTranscript(superlocusplus)
				save_winners(superlocusplus)
			superlocusplus = Superlocus(tmp_transcript)
	
	else:
		if superlocusminus == None:			# first time only
			superlocusminus = Superlocus(tmp_transcript)
		elif tmp_transcript.overlapping2(superlocusminus) >= args.fuzziness:
			superlocusminus.addTranscript(tmp_transcript)
		else:
			if superlocusminus is not None:
				superlocusminus.finishSuperlocus()
				pickBestTranscript(superlocusminus)
				save_winners(superlocusminus)
			superlocusminus = Superlocus(tmp_transcript)
		
#print("This is the last transcript processed:")
#tmp_transcript.printTranscript()
	
if superlocusplus is not None:
	superlocusplus.finishSuperlocus()
	pickBestTranscript(superlocusplus)
	save_winners(superlocusplus)
	
if superlocusminus is not None:
	superlocusminus.finishSuperlocus()
	pickBestTranscript(superlocusminus)
	save_winners(superlocusminus)

print("")
print("Done. Outputting summary (note these are rough stats, as genes can belong to multiple categories):")
stats_str = ["identical","highly similar","somewhat similar","not similar/novel"]
print("-------------------")
print("Coding genes:")
for i in range(len(stats_str)):
	print(" -",coding_stats[i],"genes are",stats_str[i],"to original annots")
	
print("-------------------")
print("Noncoding genes:")
for i in range(len(stats_str)):
	print(" -",noncoding_stats[i],"genes are",stats_str[i],"to original annots")
	
print("-------------------")
print("Times original transcripts were fused in reannotation:",merged)
print("Times original transcripts were split in reannotation:",ssplit)
print("")

print("Outputting results to output files...")

to_output = {}
def add_to_output(tt):
	global to_output
	
	if tt.chr not in to_output:
		to_output[tt.chr] = {}
	leftl = math.floor(tt.mRNA[0]/100000)
	if leftl not in to_output[tt.chr]:
		to_output[tt.chr][leftl] = []
	to_output[tt.chr][leftl].append(tt)


# output all genes while dropping antisense genes that overlap by more than 20% of either gene's length
for chr in chrlist:
#	print("Outputting annotations on pos strand of",chr)
	for loc1 in all_winners[chr+'+']:
	#	print("Loc1 is:",loc1,"corresponding to",loc1*100000,"bp")
		for fwd in all_winners[chr+'+'][loc1]:
	#		print("Looking for transcripts overlapping",fwd.finalname,"on - strand:")
	#		print(fwd.finalname,"starts at:",fwd.mRNA[0])
			if ((loc1+1)*100000) - fwd.mRNA[0] < 10000:			# transcript is close to right edge of this loc, overlapping annot may be in the next loc
				loc2rhs = loc1 + 2
			else:
				loc2rhs = loc1 + 1			
			for loc2 in range(loc1,loc2rhs):
	#			print("Looking in for antisense in loc2",loc2)
				if chr+'-' in all_winners:
					if loc2 in all_winners[chr+'-']:
						for rev in all_winners[chr+'-'][loc2]:
		#					print("Does",rev.finalname,"overlap",fwd.finalname,"?")
							overl = fwd.overlapping(rev) 
							if overl / (fwd.mRNA[1]-fwd.mRNA[0]+1) > 0.4 or overl / (rev.mRNA[1]-rev.mRNA[0]+1) > 0.4:
								# overlapping pair; which to keep?
		#						print(fwd.finalname,"and",rev.finalname,"overlap! Which to keep?")
								if fwd.coding == True and rev.coding == False:
		#							print("Keeping",fwd.finalname,"since it is coding")
									rev.is_antisense = True
								elif rev.coding == True and fwd.coding == False:
		#							print("Keeping",rev.finalname,"since it is coding")
									fwd.is_antisense = True
		#						elif fwd.coding == True and rev.coding == True:
		#							print("Both are coding, keep both")
								elif fwd.coding == False and rev.coding == False:
									if fwd.avg_cov > rev.avg_cov:
		#								print("Both are noncoding; keeping",fwd.finalname,"since it has higher expression")
										rev.is_antisense = True
									else:
		#								print("Both are noncoding; keeping",rev.finalname,"since it has higher expression")
										fwd.is_antisense = True
			add_to_output(fwd)
	
# all all remaining annots on neg strand
for chr in chrlist:
	if chr+'-' in all_winners:
		for loc1 in all_winners[chr+'-']:
			for rev in all_winners[chr+'-'][loc1]:
				add_to_output(rev)
	
			
# sort output by LHS position but no longer by strand, so output will be sorted correctly
for chr in to_output:
	for loc in to_output[chr]:
		to_output[chr][loc].sort(key = lambda x: x.mRNA[0])
		for tt in to_output[chr][loc]:
#			print("Outputting transcript:",tt.finalname)
			tt.writeToGTF(outf)
			if tt.is_antisense == False:
				tt.writeToGTF(outf_noAS)
#			else:
#				print("Is antisense, omitting from output")
		


outf.close()
outf_noAS.close()



print("Done.")
print("Full set of annotations saved to",args.outprefix+'_final.gtf')
print("Version with antisense noncoding genes removed:",args.outprefix+'_final_noAS.gtf')
print("Have a nice day!")


















