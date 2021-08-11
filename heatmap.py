

from collections import defaultdict
import HTSeq
import numpy 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from sys import argv

script, input_gtf, input_bam, output_basename = argv

minus_range = -80
plus_range = 30
largest_readlength = 98
smallest_readlength = 18

# BLOCK 1 finds valid translation start sites in the given gtf file
# criteria:
# annotated start_codon
# gene_biotype: protein_coding
# transcript_biotype: protein_coding
# transcript_support_level: 1, or 2 if no transcripts with transcript_support_level 1
# considering the above only one such start_codon per gene and only one gene and transcript per such start_codon

gtf = HTSeq.GFF_Reader(input_gtf) # HTSeq expectes 1-based coordinates in a gtf file and turns them into 0-based. Likewise intervals are turned into end-not-included.
genes_TLS = defaultdict(set)
genes_TLS2 = defaultdict(set)
TLS_genes = defaultdict(set)	
TLS_transcripts = defaultdict(list)	   
# considering genes with at least one transcript with transcript_support_level 1
for g in gtf:
    if g.type == 'start_codon' and g.attr['gene_biotype'] == 'protein_coding' and g.attr['transcript_biotype'] == 'protein_coding':
        if 'transcript_support_level' in g.attr and g.attr['transcript_support_level'] == '1':
            genes_TLS[g.attr['gene_id']].add(HTSeq.GenomicPosition(g.iv.chrom, g.iv.start_d, g.iv.strand))
# considering genes with no transcript with transcript_support_level 1, but with at least one transcript with transcript_support_level 2       
for g in gtf:
    if g.type == 'start_codon' and g.attr['gene_id'] not in genes_TLS and g.attr['gene_biotype'] == 'protein_coding' and g.attr['transcript_biotype'] == 'protein_coding':
        if 'transcript_support_level' in g.attr and g.attr['transcript_support_level'] == '2':
            genes_TLS2[g.attr['gene_id']].add(HTSeq.GenomicPosition(g.iv.chrom, g.iv.start_d, g.iv.strand))

for_del = set()
for k in genes_TLS:
    if len(genes_TLS[k]) > 1:
        for_del.add(k)
for k in for_del:
    del genes_TLS[k]  
for_del = set()
for k in genes_TLS2:
    if len(genes_TLS2[k]) > 1:
        for_del.add(k)
for k in for_del:
    del genes_TLS2[k]          
    
for g in gtf:
    if g.type == 'start_codon' and g.attr['gene_id'] in genes_TLS and g.attr['transcript_biotype'] == 'protein_coding':
        if 'transcript_support_level' in g.attr and g.attr['transcript_support_level'] == '1':
            TLS_genes[HTSeq.GenomicPosition(g.iv.chrom, g.iv.start_d, g.iv.strand)].add(g.attr['gene_id'])
            TLS_transcripts[HTSeq.GenomicPosition(g.iv.chrom, g.iv.start_d, g.iv.strand)].append(g.attr['transcript_id'])
for g in gtf:
    if g.type == 'start_codon' and g.attr['gene_id'] in genes_TLS2 and g.attr['transcript_biotype'] == 'protein_coding':
        if 'transcript_support_level' in g.attr and g.attr['transcript_support_level'] == '2':
            TLS_genes[HTSeq.GenomicPosition(g.iv.chrom, g.iv.start_d, g.iv.strand)].add(g.attr['gene_id'])
            TLS_transcripts[HTSeq.GenomicPosition(g.iv.chrom, g.iv.start_d, g.iv.strand)].append(g.attr['transcript_id'])

print len(genes_TLS)+len(genes_TLS2), 'protein_coding genes with an annotated translation start site (TLS) that is the sole annotated TLS for that gene (only considering transcripts with transcript_biotype protein_coding and transcript_support_level 1/2).'

for_del = set()
for_del2 = set()
transcript_exons = {}
for k in TLS_genes:
    if len(TLS_genes[k]) > 1:
        for_del.add(k)
    elif len(set(TLS_transcripts[k])) > 1:
        for_del2.add(k)
    else:
        transcript_exons[TLS_transcripts[k][0]]=[]
for k in for_del:
    del TLS_transcripts[k]
print 'Of that, for', len(TLS_transcripts), 'genes their TLS is unique.'
for k in for_del2:
    del TLS_transcripts[k]
print 'Of that', len(TLS_transcripts), 'TLSes with only one transcript.'

for g in gtf:
    if g.type == 'exon' and g.attr['transcript_id'] in transcript_exons:
       transcript_exons[g.attr['transcript_id']].append((g.iv.start,g.iv.end))
Starts = defaultdict(list)
for k in TLS_transcripts:
    T = TLS_transcripts[k][0]
    Starts[k] = []
    if k.strand == '+':
	    for e in range(0, len(transcript_exons[T])):
			for p in range(transcript_exons[T][e][0], transcript_exons[T][e][1]):
			    Starts[k].append(p)
    if k.strand == '-':
        for e in range(len(transcript_exons[T])-1, -1, -1):
		    for p in range(transcript_exons[T][e][1]-1, transcript_exons[T][e][0]-1, -1):
			      Starts[k].append(p)
UTR5_short = 0
CDS_short = 0				
for k in Starts:
    if Starts[k].index(k.pos) < minus_range*-1 and Starts[k].index(k.pos) != 0:
        UTR5_short += 1
    if len(Starts[k])  - Starts[k].index(k.pos) <= plus_range:
        CDS_short += 1
print UTR5_short, 'TLS/transcripts have a 5\' UTR annotated that is shorter than the desired range for plotting:', minus_range*-1, 'nt. Those are still considered.'
print CDS_short, 'TLS/transcripts have a CDS that is shorter than the desired range for plotting:', plus_range, 'nt. Those are still considered.'



# BLOCK 2 makes a GenomicArray from the input_bam		

bam = HTSeq.BAM_Reader(input_bam)
readlength_range = largest_readlength - smallest_readlength +1
coverage_readlength = HTSeq.GenomicArrayOfSets('auto', stranded=True)   
counter = 0           
for a in bam:                                                                  
 	readlength = 0
 	for b in a.cigar:
		if b.type == 'M' or b.type == 'I':
			readlength += b.size
	if readlength != 0 and largest_readlength >= readlength >= smallest_readlength:
		if coverage_readlength[a.iv.start_d_as_pos] == set([]):  
			coverage_readlength[a.iv.start_d_as_pos] = [0]*readlength_range
			coverage_readlength[a.iv.start_d_as_pos][readlength-smallest_readlength] += 1
		else:
			coverage_readlength[a.iv.start_d_as_pos][readlength-smallest_readlength] += 1
	else:
		counter += 1                                                        
print counter, 'reads had length smaller', smallest_readlength, 'or larger', largest_readlength
	
	
# BLOCK 3 collects readlength distribution for every position in the desired range for all 'approved' start sites

width = (minus_range*-1)+plus_range+1                            
t_array = numpy.zeros((width,readlength_range)) 
short_transcript = 0
for k in Starts:                     
	Chrom = k.chrom
	Pos = k.pos
	Strand = k.strand
	t_array_position = minus_range*-1  # the t_array will be filled starting from the start site going uptream and then downstream
	for p in range(0, len(Starts[k])):
		if Pos == Starts[k][p]:      
			if p >= minus_range*-1:		
				for n in range(0, (minus_range*-1)+1):  
					if coverage_readlength[HTSeq.GenomicPosition(Chrom, Starts[k][p-n], Strand)] != set([]):                
						t_array[t_array_position-n,]+= coverage_readlength[HTSeq.GenomicPosition(Chrom, Starts[k][p-n], Strand)]
			else:
			    for n in range(0, p+1):  
					if coverage_readlength[HTSeq.GenomicPosition(Chrom, Starts[k][p-n], Strand)] != set([]):                
						t_array[t_array_position-n,]+= coverage_readlength[HTSeq.GenomicPosition(Chrom, Starts[k][p-n], Strand)]
			for n in range(1, plus_range+1):  
			    if len(Starts[k]) > p+n:
				    if coverage_readlength[HTSeq.GenomicPosition(Chrom, Starts[k][p+n], Strand)] != set([]):                
					    t_array[t_array_position+n,]+= coverage_readlength[HTSeq.GenomicPosition(Chrom, Starts[k][p+n], Strand)]
			

# BLOCK 4 writes output	

numpy.savetxt(output_basename+'_start.tbl', t_array)
	
t_array_noZero = t_array.clip(0.1) # all zero counts will be converted to 0.1
t_array_log = numpy.log10(t_array_noZero)	
x=[]                                                              
for i in range(minus_range,plus_range+1):
	for r in range(len(t_array[0,])):
		x.append(i)
y=[]
for r in range(0,len(t_array)):
	for i in range(smallest_readlength,largest_readlength+1):
		y.append(i)
counts =[]
for i in range(0,len(t_array)):
	counts += list(t_array_log[i,])      

# plot the start codon heatmap
colormap = cm.jet
colormap.set_under('w')	       # color 'white' will be used for all fields below the minimum set in (plt.clim)
plt.figure(figsize=(20,10), dpi=300)
plt.hist2d(x,y, bins=[len(t_array),len(t_array[0,])], weights=counts, cmap=colormap)
plt.title('Read length distribution at positions around the annotated AUG start codon.')
plt.colorbar()
plt.clim(math.log(1,10),max(counts)) 		# the minimum is set to log10(1); all counts below 1 will be plotted as white (see line above)
plt.xticks(numpy.arange(min(x), max(x)+1, 10.0), fontsize=8, rotation=-30)
plt.ylabel('Read length [nt]')
plt.xlabel('Distance from annotated AUG start codon (A=0) [nt]') 
plt.savefig(output_basename+'_start.png')














	
