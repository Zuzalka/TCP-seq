#!/usr/bin/python

# Readlength Distribution of reads in specific areas of a transcript
# checks only gene_biotype protein-coding in the gtf
# the longest 5UTR and 3UTR of all transcripts associated with a gene is taken
# only genes are considered that have exactly one start site and exactly one stop site annotated
# if more than one gene share the same start and stop, the longest UTRs are used
# 10nt are added to both UTRs (because I see in the data that often reads start some nt before the annotated 5UTR 5'end)
# if no UTRs are annotated 100nt is used
# This script checks the FP 5' end and 3' end postition and defines their place in relation to start and stop codons and 5 and 3 UTRs. 
# If a FP falls entirely into an intron of a gene it will be counted as belonging to the gene.
# The script outputs the readlength dist. over all genes in the input for the transcript areas: start_codon, stop_codon, 5UTR, 3UTR, CDS
# definition of transcript areas: 	start_codon: the FP 5' end is at < -3 and the FP 3' end is at > +5 of A from AUG
#												5UTR: the FP 5' end is at > transcript start and the FP 3' end is at <= +5 of A from AUG
#												stop_codon: the FP 5' end is at < -6 from U of URR and the FP 3' end is at > +2 of U from URR
#												3UTR: the FP 3' end is at < transcript end and the FP 5' end is at >= -6 of U from URR
# The definitions can be changed in BLOCK 2a and 2b, first row in output should be also changed



from collections import defaultdict
from collections import Counter
import HTSeq
import numpy 
import csv
import sys
from sys import argv

script, input_gtf, input_bam, output_basename = argv     


##### BLOCK 1a #####   analysing the annotation and defining the set of suitable genes, first: for start and stop codons

gtf = HTSeq.GFF_Reader(input_gtf) 
geneIDs = set()              
not_suitable_geneIDs = set()
positions = {}
more_than_one_start = set()
more_than_one_stop = set()
for g in gtf:
	if g.type == 'gene' and g.attr['gene_biotype'] == 'protein_coding' and g.attr['gene_id'] != 'ENSG00000178971':  # this gene is CTC1 on chrom 17, it has a stag of around 300.000 reads with length of 18nt sitting in the 3UTR (at least in 2nd_total80S (the sample where I checked it)
		geneIDs.add(g.attr['gene_id'])
		positions[g.attr['gene_id']] = [0,0,0,0]
	if g.type == 'start_codon' and g.attr['gene_biotype'] == 'protein_coding' and g.attr['transcript_biotype'] == 'protein_coding' and g.attr['gene_id'] != 'ENSG00000178971':
		if positions[g.attr['gene_id']][1] == 0: 
			positions[g.attr['gene_id']][1] = g.iv.start_d
		else:
			if positions[g.attr['gene_id']][1] != g.iv.start_d:
				more_than_one_start.add(g.attr['gene_id'])
	if g.type == 'stop_codon' and g.attr['gene_biotype'] == 'protein_coding' and g.attr['transcript_biotype'] == 'protein_coding' and g.attr['gene_id'] != 'ENSG00000178971':
		if positions[g.attr['gene_id']][2] == 0: 
			positions[g.attr['gene_id']][2] = g.iv.start_d
		else:
			if positions[g.attr['gene_id']][2] != g.iv.start_d:
				more_than_one_stop.add(g.attr['gene_id'])

no_start = set()
no_stop = set()
for k in positions:
	if positions[k][1] == 0:
		no_start.add(k)
	if positions[k][2] == 0:
		no_stop.add(k)
print
print len(geneIDs), 'genes in gtf file classified as protein_coding.'
print len(more_than_one_start), 'genes with more than one annotated start_codon in transcripts with biotype protein_coding. Genes ignored.'
print len(more_than_one_stop-more_than_one_start), 'of the left genes have more than one annotated stop_codon in transcripts with biotype protein_coding. Genes ignored.'
print len(no_start.intersection(no_stop)), 'genes have no start nor stop annotated. Genes ignored.'
print len((no_start - no_stop)-(more_than_one_start.union(more_than_one_stop))), 'of the left genes have no start, but stop. Genes ignored.'
print len((no_stop - no_start)-(more_than_one_start.union(more_than_one_stop))), 'of the left genes have no stop, but start. Genes ignored.'

more_than_one_start_stop = more_than_one_start.union(more_than_one_stop)
geneIDs = geneIDs - more_than_one_start_stop
for s in more_than_one_start_stop:
	del positions[s]
no_start_stop = no_start.union(no_stop)
no_start_stop2 = no_start_stop - more_than_one_start_stop
geneIDs = geneIDs - no_start_stop2
for s in no_start_stop2:
	del positions[s]

print len(geneIDs), 'genes left.'
print



##### BLOCK 1b #####   analysing the annotation and defining the set of suitable genes, second: for UTRs 
chrom_strand = {}
for g in gtf:
	if g.attr['gene_id'] in geneIDs:
		if g.attr['gene_id'] not in chrom_strand:
			chrom_strand[g.attr['gene_id']] = (g.iv.chrom,g.iv.strand)
		if g.type == 'UTR' and g.attr['transcript_biotype'] == 'protein_coding':
			if g.iv.strand == '+':
				if g.iv.start_d < positions[g.attr['gene_id']][1]:   # than it is 5UTR
					if positions[g.attr['gene_id']][0] == 0:
						positions[g.attr['gene_id']][0] = g.iv.start_d
					else:
						if g.iv.start_d < positions[g.attr['gene_id']][0]:
							positions[g.attr['gene_id']][0] = g.iv.start_d
				if g.iv.start_d > positions[g.attr['gene_id']][2]:   # than it is 3UTR
					if positions[g.attr['gene_id']][3] == 0:
						positions[g.attr['gene_id']][3] = g.iv.end_d-1
					else:
						if g.iv.end_d-1 > positions[g.attr['gene_id']][3]:
							positions[g.attr['gene_id']][3] = g.iv.end_d-1
			if g.iv.strand == '-':
				if g.iv.start_d > positions[g.attr['gene_id']][1]:   # than it is 5UTR
					if positions[g.attr['gene_id']][0] == 0:
						positions[g.attr['gene_id']][0] = g.iv.start_d
					else:
						if g.iv.start_d > positions[g.attr['gene_id']][0]:
							positions[g.attr['gene_id']][0] = g.iv.start_d
				if g.iv.start_d < positions[g.attr['gene_id']][2]:   # than it is 3UTR
					if positions[g.attr['gene_id']][3] == 0:
						positions[g.attr['gene_id']][3] = g.iv.end_d+1
					else:
						if g.iv.end_d+1 < positions[g.attr['gene_id']][3]:
							positions[g.attr['gene_id']][3] = g.iv.end_d+1



#### BLOCK 1c #### analysing the annotation and defining the set of suitable genes, third: remove genes with same start/stop but different stop/start; for genes with same start and stop choose longest 5UTR and 3UTR

starts = []
stops = []
for k in positions:
	starts.append(chrom_strand[k][0]+'/'+chrom_strand[k][1]+'/'+str(positions[k][1]))
	stops.append(chrom_strand[k][0]+'/'+chrom_strand[k][1]+'/'+str(positions[k][2]))
count_starts = Counter(starts)
count_stops = Counter(stops)
starts = []
stops = []
for c in count_starts:
	starts.append(count_starts[c])
for c in count_stops:
	stops.append(count_stops[c])
print 'That often a start is shared by that many genes:'
print Counter(starts)
print 'That often a stop is shared by that many genes:'
print Counter(stops)
print

count_same_start_different_stop = 0
for c in count_starts:
	if count_starts[c] > 1:
		intermediate_pos = set()
		intermediate_geneIDs = []
		for k in positions:
			if chrom_strand[k][0]+'/'+chrom_strand[k][1]+'/'+str(positions[k][1]) == c:
				intermediate_pos.add(positions[k][2])
				intermediate_geneIDs.append(k)
		if len(intermediate_geneIDs) != count_starts[c]:
			print 'error'
		if len(intermediate_pos) > 1:
			count_same_start_different_stop += 1
			for k in intermediate_geneIDs:
				del positions[k]
				geneIDs.remove(k)
		if len(intermediate_pos) == 1:
			UTR5 = []
			UTR3 = []
			for k in intermediate_geneIDs:
				if positions[k][0] != 0:
					UTR5.append(positions[k][0])
				if positions[k][3] != 0:
					UTR3.append(positions[k][3])
			if c.split('/')[1] == '+':
				if len(UTR5) != 0:
					positions[intermediate_geneIDs[0]][0] = min(UTR5)
				if len(UTR3) != 0:
					positions[intermediate_geneIDs[0]][3] = max(UTR3)
			if c.split('/')[1] == '-':
				if len(UTR5) != 0:
					positions[intermediate_geneIDs[0]][0] = max(UTR5)
				if len(UTR3) != 0:
					positions[intermediate_geneIDs[0]][3] = min(UTR3)
			for i in range(1, len(intermediate_geneIDs)):
				del positions[intermediate_geneIDs[i]]   # only one entry is kept of all genes that have the same stop and start, the longest 5' and 3'UTR is choosen, or if none annotated the value stays 0
				
print count_same_start_different_stop, ' times the genes that share a start have not all the same stop. Genes ignored.'				
print len(positions), 'genes left.'

count_same_stop_different_start = 0
for c in count_stops:
	if count_stops[c] > 1:
		intermediate_pos = set()
		intermediate_geneIDs = []
		for k in positions:
			if chrom_strand[k][0]+'/'+chrom_strand[k][1]+'/'+str(positions[k][2]) == c:
				intermediate_pos.add(positions[k][1])
				intermediate_geneIDs.append(k)
		if len(intermediate_pos) > 1:                    # All genes with same start and stop were found in the previous module. Here all genes are removed that have the same stop but different start (therefore len(intermediate_pos) can't be == 1)
			count_same_stop_different_start += 1
			for k in intermediate_geneIDs:
				del positions[k]
				geneIDs.remove(k)
						
print count_same_stop_different_start, ' times the genes (of the left genes) that share a stop have not all the same start. Genes ignored.'				
print len(positions), 'genes left.'



#### BLOCK 2a ####	    defining stop and start (and therefore UTR and CDS) intervals taking introns in account

coverage_all_exons = HTSeq.GenomicArray('auto', stranded=True, typecode='i') 
for g in gtf:
	if g.attr['gene_id'] in geneIDs:     # in geneIDs are all considered genes (also all geneIDs that share a start and stop)
		if g.type == 'exon' and g.attr['transcript_biotype'] == 'protein_coding':	
			coverage_all_exons[g.iv] += 1

i_positions = {}
upsZero = 0
for k in positions:   
	if chrom_strand[k][1] == '+':
		istart_minus = 0
		istart_plus = 0
		if positions[k][0] == 0 or positions[k][1] - positions[k][0] < 3:   # that means if no or too little 5UTR is annotated
			istart_minus = positions[k][1]-3
		else:	
			n = 1
			exon_dist = 0
			while True:
				if coverage_all_exons[HTSeq.GenomicPosition(chrom_strand[k][0], positions[k][1]-n, chrom_strand[k][1])] != 0:
					n += 1
					exon_dist -= 1
				else:
					n += 1
				if exon_dist == -3:
					istart_minus = positions[k][1]-(n-1)
					break 
					
		n = 1
		exon_dist = 0
		while True:
			if coverage_all_exons[HTSeq.GenomicPosition(chrom_strand[k][0], positions[k][1]+n, chrom_strand[k][1])] != 0:
				n += 1
				exon_dist += 1
			else:
				n += 1
			if exon_dist == 5:
				istart_plus = positions[k][1]+(n-1)
				break
		
		istop_minus = 0
		istop_plus = 0 
		n = 1
		exon_dist = 0
		while True:
			if coverage_all_exons[HTSeq.GenomicPosition(chrom_strand[k][0], positions[k][2]-n, chrom_strand[k][1])] != 0:
				n += 1
				exon_dist -= 1
			else:
				n += 1
			if exon_dist == -6:
				istop_minus = positions[k][2]-(n-1)
				break 
		
#		if positions[k][3] == 0 or positions[k][3] - positions[k][2] < 2:   # that means if no or too little 3UTR is annotated
#			istop_plus = positions[k][2]+2
#		else:				
		n = 1
		exon_dist = 0
		while True:
			if coverage_all_exons[HTSeq.GenomicPosition(chrom_strand[k][0], positions[k][2]+n, chrom_strand[k][1])] != 0:
				n += 1
				exon_dist += 1
			else:
				n += 1
			if exon_dist == 2:
				istop_plus = positions[k][2]+(n-1)
				break
						
		if istart_minus != 0 and istart_plus != 0 and istop_minus != 0 and istop_plus != 0:
			i_positions[k] = [positions[k][0],istart_minus,istart_plus,istop_minus,istop_plus,positions[k][3]]
		else:
			upsZero += 1
		
	if chrom_strand[k][1] == '-':
		istart_minus = 0
		istart_plus = 0 
		if positions[k][0] == 0 or positions[k][0] - positions[k][1] < 3:   # that means if no or too little 5UTR is annotated
			istart_minus = positions[k][1]+3
		else:	
			n = 1
			exon_dist = 0
			while True:
				if coverage_all_exons[HTSeq.GenomicPosition(chrom_strand[k][0], positions[k][1]+n, chrom_strand[k][1])] != 0:
					n += 1
					exon_dist -= 1
				else:
					n += 1
				if exon_dist == -3:
					istart_minus = positions[k][1]+(n-1)
					break 
			
		n = 1
		exon_dist = 0
		while True:
			if coverage_all_exons[HTSeq.GenomicPosition(chrom_strand[k][0], positions[k][1]-n, chrom_strand[k][1])] != 0:
				n += 1
				exon_dist += 1
			else:
				n += 1
			if exon_dist == 5:
				istart_plus = positions[k][1]-(n-1)
				break
		
		istop_minus = 0
		istop_plus = 0 		
		n = 1
		exon_dist = 0
		while True:
			if coverage_all_exons[HTSeq.GenomicPosition(chrom_strand[k][0], positions[k][2]+n, chrom_strand[k][1])] != 0:
				n += 1
				exon_dist -= 1
			else:
				n += 1
			if exon_dist == -6:
				istop_minus = positions[k][2]+(n-1)
				break 
				
#		if positions[k][3] == 0 or positions[k][2] - positions[k][3] < 2:   # that means if no 3UTR is annotated
#			istop_plus = positions[k][2]-2
#		else:	
		n = 1
		exon_dist = 0
		while True:
			if coverage_all_exons[HTSeq.GenomicPosition(chrom_strand[k][0], positions[k][2]-n, chrom_strand[k][1])] != 0:
				n += 1
				exon_dist += 1
			else:
				n += 1
			if exon_dist == 2:
				istop_plus = positions[k][2]-(n-1)
				break		
		
		if istart_minus != 0 and istart_plus != 0 and istop_minus != 0 and istop_plus != 0:
			i_positions[k] = [positions[k][0],istart_minus,istart_plus,istop_minus,istop_plus,positions[k][3]]
		else:
			upsZero += 1		
					
print upsZero, 'This number should be zero.'


#### BLOCK 2b ####    add 100nt UTRs if no UTRs are annotated and 10nt to annotated UTRs
count_no_5and3UTRs = 0
count_no_5UTR = 0
count_no_3UTR = 0
for k in i_positions:
	if chrom_strand[k][1] == '+':
		if i_positions[k][0] == 0 and i_positions[k][5] == 0:
			count_no_5and3UTRs += 1
			i_positions[k][0] = i_positions[k][1] - (100-3)
			i_positions[k][5] = i_positions[k][4] + (100-2)
		if i_positions[k][0] == 0 and i_positions[k][5] != 0:
			count_no_5UTR += 1
			i_positions[k][0] = i_positions[k][1] - (100-3)
			i_positions[k][5] += 10
		if i_positions[k][0] != 0 and i_positions[k][5] == 0:
			count_no_3UTR += 1
			i_positions[k][0] -= 10
			i_positions[k][5] = i_positions[k][4] + (100-2)
	if chrom_strand[k][1] == '-':
		if i_positions[k][0] == 0 and i_positions[k][5] == 0:
			count_no_5and3UTRs += 1
			i_positions[k][0] = i_positions[k][1] + (100-3)
			i_positions[k][5] = i_positions[k][4] - (100-2)
		if i_positions[k][0] == 0 and i_positions[k][5] != 0:
			count_no_5UTR += 1
			i_positions[k][0] = i_positions[k][1] + (100-3)
			i_positions[k][5] -= 10
		if i_positions[k][0] != 0 and i_positions[k][5] == 0:
			count_no_3UTR += 1
			i_positions[k][0] += 10
			i_positions[k][5] = i_positions[k][4] - (100-2)
print
print count_no_5and3UTRs, ' of the left genes have no 5\' nor 3\' UTR annotated (100nt used).'
print count_no_5UTR, ' of the left genes have no 5\'UTR annotated (100nt used).'
print count_no_3UTR, ' of the left genes have no 3\'UTR annotated (100nt used).'
print 'To annotated UTRs 10nt were added.'


### BLOCK 3 ####   changing data structure of i_positions

final = defaultdict(list)
for k in i_positions:
	final[chrom_strand[k][0]+chrom_strand[k][1]].append(i_positions[k])


#### BLOCK 4 ####    readlength distribution

bam = HTSeq.BAM_Reader(input_bam)
smallest_readlength = 18 
largest_readlength = 100 
readlength_range = largest_readlength - smallest_readlength +1
RLdist_5UTR = [0] * readlength_range	
RLdist_3UTR = [0] * readlength_range
RLdist_start = [0] * readlength_range
RLdist_stop = [0] * readlength_range
RLdist_CDS = [0] * readlength_range

no_place = 0
one_place = 0
no_of_places = []  
cigar_types = set() 
more_places_in_gene = [] 
read_in_intron = 0  
for a in bam:
	counter=0
	readlength = 0
	for b in a.cigar:
		cigar_types.add(b.type)
		if b.type == 'M' or b.type == 'I':  # cigars N and D have no influence on the readlength
			readlength += b.size
	if a.iv.strand == '+' and 18 <= readlength <= 100:
		if a.iv.chrom+a.iv.strand in final:
			for k in final[a.iv.chrom+a.iv.strand]:   
	 			if a.iv.start_d >= k[0] and a.iv.end_d-1 <= k[5]:     	# read belongs to a gene
	 				counter+=1
	 				counter2 = 0
	 				counter_intron = 0
	 				if a.iv.end_d-1 <= k[2]: 										# belongs to 5UTR
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:				
	 						RLdist_5UTR[readlength-smallest_readlength]+=1
	 						counter2 += 1
	 					else: 
	 						counter_intron += 1
	 				if a.iv.start_d < k[1] and a.iv.end_d-1 > k[2]:     	# belongs to start
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:	
	 						RLdist_start[readlength-smallest_readlength]+=1
	 						counter2 += 1
	 					else: 
	 						counter_intron += 1	
	 				if a.iv.start_d >= k[1] and a.iv.end_d-1 <= k[4]:    	# belongs to CDS
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:	
	 						RLdist_CDS[readlength-smallest_readlength]+=1
	 						counter2 += 1
	 					else: 
	 						counter_intron += 1	
	 				if a.iv.start_d < k[3] and a.iv.end_d-1 > k[4]:     	# belongs to stop
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:	
	 						RLdist_stop[readlength-smallest_readlength]+=1
	 						counter2 += 1	
	 					else: 
	 						counter_intron += 1	
	 				if a.iv.start_d >= k[3]:										# belongs to 3UTR
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:	
	 						RLdist_3UTR[readlength-smallest_readlength]+=1
	 						counter2 += 1
	 					else: 
	 						counter_intron += 1	
	 				more_places_in_gene.append(counter2)
					if counter2 == 0 and counter_intron > 0:
						read_in_intron += 1

	
	if a.iv.strand == '-' and 18 <= readlength <= 100:
		if a.iv.chrom+a.iv.strand in final:
			for k in final[a.iv.chrom+a.iv.strand]:   
	 			if a.iv.start_d <= k[0] and a.iv.end_d+1 >= k[5]:     	# read belongs to a gene
	 				counter+=1
	 				counter2 = 0
	 				counter_intron = 0
	 				if a.iv.end_d+1 >= k[2]:     								 	# belongs to 5UTR
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:	
	 						RLdist_5UTR[readlength-smallest_readlength]+=1
	 						counter2 += 1
	 					else: 
	 						counter_intron += 1		
	 				if a.iv.start_d > k[1] and a.iv.end_d+1 < k[2]:     	# belongs to start
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:	
	 						RLdist_start[readlength-smallest_readlength]+=1
	 						counter2 += 1
	 					else: 
	 						counter_intron += 1		
	 				if a.iv.start_d <= k[1] and a.iv.end_d+1 >= k[4]:    	# belongs to CDS
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:	
	 						RLdist_CDS[readlength-smallest_readlength]+=1
	 						counter2 += 1
	 					else: 
	 						counter_intron += 1		
	 				if a.iv.start_d > k[3] and a.iv.end_d+1 < k[4]:     	# belongs to stop
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:	
	 						RLdist_stop[readlength-smallest_readlength]+=1
	 						counter2 += 1
	 					else: 
	 						counter_intron += 1			
	 				if a.iv.start_d <= k[3]:										# belongs to 3UTR
	 					count_nt = 0
	 					count_exon = 0
	 					for b in a.cigar:  
	 						if b.type == 'M':
	 							for nt in b.ref_iv.xrange(step=1):
	 								count_nt +=1 
	 								if coverage_all_exons[nt] != 0:
	 									count_exon +=1
	 					if count_exon == count_nt:	
	 						RLdist_3UTR[readlength-smallest_readlength]+=1
	 						counter2 += 1
	 					else: 
	 						counter_intron += 1		
	 				more_places_in_gene.append(counter2)
					if counter2 == 0 and counter_intron > 0:
						read_in_intron += 1
	 					
	if counter == 0:
	 	no_place += 1
	if counter == 1:
	 	one_place += 1
	if counter >1:
	 	no_of_places.append(counter)


d = Counter(more_places_in_gene)

print 'That often an alignment found that many places in one gene. It should be only one.'
print d
print read_in_intron, 'times an alignment was placed in an intron of a gene. Not counted.'	
	
	
c=Counter(no_of_places)

print
print cigar_types, 'cigar types in bam_input. (M and I used for read lenght calculation)'
print                                                   
print one_place, 'times an alignment was counted for exactly one gene.'
print no_place, 'time an alignment was counted for none of the genes.'
print len(no_of_places), 'times an alignment was counted for more than one gene. (ideally, this should be zero) (all places are counted)'
print 'list of how often an alignment was counted for that many genes:'
print 'no. of genes         for how many alignments'
for i in c:
	print i, c[i]
	
	
	
	
#### BLOCK 5 ####      output
output_filename = output_basename+'_RLdist.tbl'
with open(output_filename, 'wb') as o:
		out = csv.writer(o, delimiter='\t')
		out.writerow([script, 'start interval from -3 to +5 of A of AUG', 'stop interval from -6 to +2 of U of URR'] )
		out.writerow(['readlength', 'count 5UTR','count start','count CDS','count stop','count 3UTR'])
		for j in range(0, readlength_range):
			out.writerow([j+smallest_readlength,RLdist_5UTR[j],RLdist_start[j],RLdist_CDS[j],RLdist_stop[j],RLdist_3UTR[j]])	
		out.writerow(['sum',sum(RLdist_5UTR),sum(RLdist_start),sum(RLdist_CDS),sum(RLdist_stop),sum(RLdist_3UTR)])

	




