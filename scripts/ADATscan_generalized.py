#!/usr/bin/env python3
"""
python ADATscan_generalized.py <fasta> <null model output> <Enrichment and Depletion output> <gene body output>

user input the codons of interest 
 
"""


#Import modules needed 
import sys
import numpy as np
import scipy
from scipy import stats

codon_set = []


while 1 == 1:
	x = input('Input codon of interest, if complete then type: \'n\'')
	if x.lower() == 'n': break
	if x.lower() == '\'n\'': break
	if len(x) != 3: 
		print('Error: Improper codon length, must be three bases!')
		continue
	if len(x) == 3 and x.count('A') + x.count('a') + x.count('T') + x.count('t') + x.count('C') + x.count('c') + x.count('G') + x.count('g') != 3: 
		print ('Error: Improper codon composition, must be A/T/C/G only!')
		continue
	#reject ATG and TGG and TGA and TAG and TAA because there cannot be bias for those in this way 
	if x.upper() == 'ATG': 
		print('Invalid codon, methionine is only encoded by one codon')
		continue

	if x.upper() == 'TGG':
		print('Invalid codon, tryptophan is only encoded by one codon')
		continue
	
	if x.upper() == 'TGA': 
		print('Invalid codon, stop codons cannot be enriched or depleted')
		continue
	if x.upper() == 'TAG': 
		print('Invalid codon, stop codons cannot be enriched or depleted')
		continue
	if x.upper() == 'TAA': 
		print('Invalid codon, stop codons cannot be enriched or depleted')
		continue
	

	codon_set.append(x.upper())

while 1 == 1: 
	x = input("Which multiple test correction will be run? (Bonferroni (B) or Benjamani-hochberg (BH))")
	if x.upper() == 'B': 
		correction = 'bonferroni'
		break
	if x.upper() == 'BH':
		correction = 'benjamani-hochberg'
		x = input("What FDR will be used? (ex. 0.1)")
		FDR = float(x)
		break
	print('Problem with input. Type either \'B\' or \'BH\'')

		
##################STEP 1 GENERATE DATA STRUCTURES############################## 
#specify genetic code
codons = {}
codons[ "A" ] = ["GCA", "GCC", "GCG", "GCT" ] 
codons[ "C" ] = ["TGC", "TGT" ]
codons[ "D" ] = ["GAC", "GAT" ]
codons[ "E" ] = ["GAA", "GAG" ]
codons[ "F" ] = ["TTC", "TTT" ]
codons[ "G" ] = ["GGA", "GGC", "GGG", "GGT" ]
codons[ "H" ] = ["CAC", "CAT" ]
codons[ "I" ] = ["ATA", "ATC", "ATT" ]
codons[ "K" ] = ["AAA", "AAG" ]
codons[ "L" ] = ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG" ]
codons[ "M" ] = ["ATG" ]
codons[ "N" ] = ["AAC", "AAT" ]
codons[ "P" ] = ["CCA", "CCC", "CCG", "CCT" ]
codons[ "Q" ] = ["CAA", "CAG" ]
codons[ "R" ] = ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT" ]
codons[ "S" ] = ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT" ]
codons[ "T" ] = ["ACA", "ACC", "ACG", "ACT" ]
codons[ "V" ] = ["GTA", "GTC", "GTG", "GTT" ]
codons[ "W" ] = ["TGG" ]
codons[ "Y" ] = ["TAC", "TAT" ]
codons[ "*" ] = ["TAA", "TAG", "TGA" ]

#Reverse the dictionary so that codons are keys and AAs are values 
translatedict = {}
listofcodons = []
listofAAs = []
for k,v in codons.items():
 for i in range(len(v)):
  translatedict[v[i]] = [k]


aa_set = []
for i in range(len(codon_set)):
	x = translatedict[codon_set[i]][0]
	y = x in aa_set
	if y == True: continue
	aa_set.append(translatedict[codon_set[i]][0])


print(aa_set)
 #This set follows Percudanietal1997
##################STEP 1 GENERATE DATA STRUCTURES############################## 





##################STEP 2 GENERATE NULL MODEL OF CODON BIAS IN SAMPLE############################## 
#read in the protein coding seqs
hexome = open(sys.argv[1])

genes = [] #hold gene names 

#initialize variables to hold counts of aa and adat codon usage (TAPSLIVR)
tcount = 0
tadatreqcount = 0
acount = 0
aadatreqcount = 0
pcount = 0
padatreqcount = 0
scount = 0
sadatreqcount = 0
lcount = 0
ladatreqcount = 0
icount = 0
iadatreqcount = 0
vcount = 0
vadatreqcount = 0
rcount = 0
radatreqcount = 0

cadatreqcount = 0
ccount = 0
dadatreqcount = 0
dcount = 0
eadatreqcount = 0
ecount = 0
fadatreqcount = 0
fcount = 0
gadatreqcount = 0
gcount = 0
hadatreqcount = 0
hcount = 0
kadatreqcount = 0
kcount = 0
nadatreqcount = 0
ncount = 0
qadatreqcount = 0
qcount = 0
yadatreqcount = 0
ycount = 0


counter = 0 #initialize counter to count total lines in file 



start_count = 0
stop_count = 0
skipped = False 
for line in hexome:
	counter += 1 	
	line_list = line.strip().split('|') #parse the line by the character following the gene name 
	duplicate = line_list[0][1:] in genes
	if duplicate == True: continue

	if line_list[0][0] == '>': # if the line is a header line record the gene name 
		genes.append(line_list[0][1:])#add the gene name to the genes list 
		continue 
	if line_list[0][0:3] != 'ATG': #if the gene does not start with a start codon remove it from the list and continue
		genes = genes[0:(len(genes)-1)] #if wrong then truncate the list and move on 
		start_count +=1 
		continue
	x = line_list[0][len(line_list[0])-3:] 
	stop_codon = 0
	if x == 'TAA': stop_codon += 1 
	if x == 'TGA': stop_codon += 1 
	if x == 'TAG': stop_codon += 1 
	if stop_codon == 0: #if the gene does not end with a stop codon then remove it from the list and count it as filterred out 
		stop_count += 1 
		genes = genes[0:(len(genes)-1)]#if wrong then truncate the list and move on 
		continue
	for i in range(int((len(line_list[0]))/3)): #iterate through discreet windows of three nucleotides 
		current_codon = line_list[0][i*3:(i*3)+3] #grab out each codon
		detect_N = current_codon in translatedict
		if detect_N == False: continue 
		currentaa = translatedict[current_codon] #translate the codon using dictionary made above 
		
		if currentaa == ['C'] and 'C' in aa_set: ccount += 1 
		if current_codon == 'TGC' and 'TGC' in codon_set  : cadatreqcount +=1 
		if current_codon == 'TGT' and 'TGT' in codon_set  : cadatreqcount +=1 
		if currentaa == ['D'] and 'D' in aa_set: dcount += 1 
		if current_codon == 'GAC' and 'GAC' in codon_set  : dadatreqcount +=1 
		if current_codon == 'GAT' and 'GAT' in codon_set  : dadatreqcount +=1 
		if currentaa == ['E'] and 'E' in aa_set: ecount += 1 
		if current_codon == 'GAA' and 'GAA' in codon_set  : eadatreqcount +=1 
		if current_codon == 'GAG' and 'GAG' in codon_set  : eadatreqcount +=1 
		if currentaa == ['F'] and 'F' in aa_set: fcount += 1 
		if current_codon == 'TTC' and 'TTC' in codon_set  : fadatreqcount +=1 
		if current_codon == 'TTT' and 'TTT' in codon_set  : fadatreqcount +=1 
		if currentaa == ['G'] and 'G' in aa_set: gcount += 1 
		if current_codon == 'GGA' and 'GGA' in codon_set  : gadatreqcount +=1 
		if current_codon == 'GGC' and 'GGC' in codon_set  : gadatreqcount +=1 
		if current_codon == 'GGG' and 'GGG' in codon_set  : gadatreqcount +=1 
		if current_codon == 'GGT' and 'GGT' in codon_set  : gadatreqcount +=1 
		if currentaa == ['H'] and 'H' in aa_set: hcount += 1 
		if current_codon == 'CAC' and 'CAC' in codon_set  : hadatreqcount +=1 
		if current_codon == 'CAT' and 'CAT' in codon_set  : hadatreqcount +=1 
		if currentaa == ['Y'] and 'Y' in aa_set: hcount += 1 
		if current_codon == 'TAC' and 'TAC' in codon_set  : yadatreqcount +=1 
		if current_codon == 'TAT' and 'TAT' in codon_set  : yadatreqcount +=1 
		if currentaa == ['K'] and 'K' in aa_set: kcount += 1 
		if current_codon == 'AAA' and 'AAA' in codon_set  : kadatreqcount +=1 
		if current_codon == 'AAG' and 'AAG' in codon_set : kadatreqcount +=1 
		if currentaa == ['N'] and 'N' in aa_set: ncount += 1 
		if current_codon == 'AAC' and 'AAC' in codon_set  : nadatreqcount +=1 
		if current_codon == 'AAT' and 'AAT' in codon_set  : nadatreqcount +=1 
		if currentaa == ['Q'] and 'Q' in aa_set: qcount += 1 
		if current_codon == 'CAA' and 'CAA' in codon_set  : qadatreqcount +=1 
		if current_codon == 'CAG' and 'CAG' in codon_set : qadatreqcount +=1 
		if currentaa == ['T'] and 'T' in aa_set: tcount += 1 
		if current_codon == 'ACC' and 'ACC' in codon_set  : tadatreqcount +=1 
		if current_codon == 'ACA' and 'ACA' in codon_set  : tadatreqcount +=1 
		if current_codon == 'ACG' and 'ACG' in codon_set : tadatreqcount +=1 
		if current_codon == 'ACT' and 'ACT' in codon_set  : tadatreqcount +=1 
		if currentaa == ['A'] and 'A' in aa_set: acount += 1 
		if current_codon == 'GCC' and 'GCC' in codon_set : aadatreqcount +=1 
		if current_codon == 'GCA' and 'GCA' in codon_set : aadatreqcount +=1 
		if current_codon == 'GCG' and 'GCG' in codon_set : aadatreqcount +=1 
		if current_codon == 'GCT' and 'GCT' in codon_set : aadatreqcount +=1 
		if currentaa == ['P'] and 'P' in aa_set: pcount += 1 
		if current_codon == 'CCC' and 'CCC' in codon_set : padatreqcount +=1 
		if current_codon == 'CCA' and 'CCA' in codon_set : padatreqcount +=1 
		if current_codon == 'CCG' and 'CCG' in codon_set : padatreqcount +=1 
		if current_codon == 'CCT' and 'CCT' in codon_set : padatreqcount +=1 
		if currentaa == ['S'] and 'S' in aa_set: scount += 1 
		if current_codon == 'TCC' and 'TCC' in codon_set: sadatreqcount +=1 
		if current_codon == 'TCA' and 'TCA' in codon_set: sadatreqcount +=1 
		if current_codon == 'TCG' and 'TCG' in codon_set: sadatreqcount +=1 
		if current_codon == 'TCT' and 'TCT' in codon_set: sadatreqcount +=1 
		if current_codon == 'AGC' and 'AGC' in codon_set: sadatreqcount +=1 
		if current_codon == 'AGT' and 'AGT' in codon_set: sadatreqcount +=1 
		if currentaa == ['L'] and 'L' in aa_set: lcount += 1 
		if current_codon == 'CTC' and 'CTC' in codon_set: ladatreqcount +=1 
		if current_codon == 'CTA' and 'CTA' in codon_set: ladatreqcount +=1 
		if current_codon == 'CTG' and 'CTG' in codon_set: ladatreqcount +=1 
		if current_codon == 'CTT' and 'CTT' in codon_set: ladatreqcount +=1 
		if current_codon == 'TTA' and 'TTA' in codon_set: ladatreqcount +=1 
		if current_codon == 'TTG' and 'TTG' in codon_set: ladatreqcount +=1 
		if currentaa == ['I'] and 'I' in aa_set: icount += 1 
		if current_codon == 'ATC' and 'ATC' in codon_set : iadatreqcount +=1 
		if current_codon == 'ATA' and 'ATA' in codon_set : iadatreqcount +=1 
		if current_codon == 'ATT' and 'ATT' in codon_set : iadatreqcount +=1 
		if currentaa == ['V'] and 'V' in aa_set: vcount += 1 
		if current_codon == 'GTC' and 'GTC' in codon_set : vadatreqcount +=1 
		if current_codon == 'GTA' and 'GTA' in codon_set : vadatreqcount +=1 
		if current_codon == 'GTG' and 'GTG' in codon_set : vadatreqcount +=1 
		if current_codon == 'GTT' and 'GTT' in codon_set : vadatreqcount +=1 
		if currentaa == ['R'] and 'R' in aa_set: rcount += 1 
		if current_codon == 'AGA' and 'AGA' in codon_set: radatreqcount +=1 
		if current_codon == 'AGG' and 'AGG' in codon_set: radatreqcount +=1
		if current_codon == 'CGC' and 'CGC' in codon_set: radatreqcount +=1 
		if current_codon == 'CGA' and 'CGA' in codon_set: radatreqcount +=1
		if current_codon == 'CGG' and 'CGG' in codon_set: radatreqcount +=1 
		if current_codon == 'CGT' and 'CGT' in codon_set: radatreqcount +=1










hexome.close()#close exome file 

#calculate genomic background of codon usage 
if 'T' in aa_set: tfrac = tadatreqcount/tcount 
else: tfrac = 0 
if 'A' in aa_set: afrac = aadatreqcount/acount
else: afrac = 0 
if 'P' in aa_set: pfrac = padatreqcount/pcount
else: pfrac = 0 
if 'S' in aa_set: sfrac = sadatreqcount/scount
else: sfrac = 0 
if 'L' in aa_set: lfrac = ladatreqcount/lcount
else: lfrac = 0 
if 'I' in aa_set: ifrac = iadatreqcount/icount
else: ifrac = 0 
if 'V' in aa_set: vfrac = vadatreqcount/vcount
else: vfrac = 0 
if 'R' in aa_set: rfrac = radatreqcount/rcount
else: rfrac = 0 
if 'C' in aa_set: cfrac = cadatreqcount/ccount
else: cfrac = 0 
if 'D' in aa_set: dfrac = dadatreqcount/dcount
else: dfrac = 0 
if 'E' in aa_set: efrac = eadatreqcount/ecount
else: efrac = 0 
if 'F' in aa_set: ffrac = fadatreqcount/fcount
else: ffrac = 0 
if 'G' in aa_set: gfrac = gadatreqcount/gcount
else: gfrac = 0 
if 'H' in aa_set: hfrac = hadatreqcount/hcount
else: hfrac = 0 	
if 'K' in aa_set: kfrac = kadatreqcount/kcount
else: kfrac = 0 
if 'N' in aa_set: nfrac = nadatreqcount/ncount
else: nfrac = 0 	
if 'Q' in aa_set: qfrac = qadatreqcount/qcount
else: qfrac = 0 
if 'Y' in aa_set: yfrac = yadatreqcount/ycount
else: yfrac = 0 






#Print general results of this upstream loop to the null model file 
original_stdout = sys.stdout 
f = open(sys.argv[2], 'w') 
sys.stdout = f
print("Genes lacking start codon count:", start_count,"Genes lacking stop codon count:", stop_count)
print("Genes retained from file:", len(genes))
print("Total genes in file:",int(counter/2))
if 'T' in aa_set: print('Codon of interest fraction threonine:'+str(tfrac))
if 'A' in aa_set: print('Codon of interest fraction alanine:'+str(afrac))
if 'P' in aa_set: print('Codon of interest fraction proline:'+str(pfrac)) 
if 'S' in aa_set: print('Codon of interest fraction serine:'+str(sfrac ))
if 'L' in aa_set: print('Codon of interest fraction luecine:'+str(lfrac ))
if 'I' in aa_set: print('Codon of interest fraction isoluecine:'+str(ifrac ))
if 'V' in aa_set: print('Codon of interest fraction valine:'+str(vfrac ))
if 'R' in aa_set: print('Codon of interest fraction arginine:'+str(rfrac)) 
if 'C' in aa_set: print('Codon of interest fraction cysteine:'+str(cfrac))
if 'D' in aa_set: print('Codon of interest fraction aspartic acid:'+str(dfrac))
if 'E' in aa_set: print('Codon of interest fraction glutamic acid:'+str(efrac)) 
if 'F' in aa_set: print('Codon of interest fraction phenylalanine:'+str(ffrac ))
if 'G' in aa_set: print('Codon of interest fraction glycine:'+str(gfrac ))
if 'H' in aa_set: print('Codon of interest fraction histidine:'+str(hfrac ))
if 'K' in aa_set: print('Codon of interest fraction lysine:'+str(kfrac ))
if 'N' in aa_set: print('Codon of interest fraction asparagine:'+str(ifrac ))
if 'Q' in aa_set: print('Codon of interest fraction glutamine:'+str(qfrac ))
if 'Y' in aa_set: print('Codon of interest fraction tyrosine:'+str(yfrac)) 



zzz = len(genes)
##################STEP 2 GENERATE NULL MODEL OF CODON BIAS IN SAMPLE############################## 





####################################################################################################################
##################STEP 3 iterate through fasta get number of ADAT-dependent codons##################################
####################################################################################################################

#read in the protein coding seqs
hexome = open(sys.argv[1])
#initialize variables to be iterated 
genes = [] #holds gene names 
filter_counter = 0 #counts number of genes filterred out 
ADAT_Codon_counts = [] #count number of ADAT codons in each gene 
ADAT_Codon_freqs = [] #record frequency of ADAT codons in each gene (normalize by length)
ADAT_pvalues_t = [] #retain p values for chisquare tests of enrichment for adat codons (transform so they can be ordered by enrichment) 
ADAT_pvalues = []
adat_exp = []
adat_obs = []
protein_length_AA = []
protein_length_nucleotides = [] #double check to ensure the number of nucleotides is right for each gene
bonferroni_corrected_chisquare = []
bonferroni_enriched_depleted = []
skipped = False

start_count = 0
stop_count = 0
gene_body_dict = {}

for line in hexome: #for line in fasta file 
	line_list = line.strip().split('|') #parse the line by the character following the gene name 
	duplicate = line_list[0][1:] in genes
	if duplicate == True: continue
	if line_list[0][0] == '>': # if the line is a header line record the gene name 
		genes.append(line_list[0][1:]) #append the gene name 
		ADAT_Codon_counts.append(0) #append a zero to be incremented later 
		protein_length_nucleotides.append(0) #append a zero to be incremented later 
		continue #move on to sequence line 
	if line_list[0][0:3] != 'ATG': #if the gene does not start with a start codon remove it from the list and count it as filterred out 
		start_count +=1 
		genes = genes[0:(len(genes)-1)] #correct way to truncate a list (remove the last item via slicing if filter not met)
		ADAT_Codon_counts = ADAT_Codon_counts[0:len(ADAT_Codon_counts)-1] #same for the variables counting stuff 
		protein_length_nucleotides = protein_length_nucleotides[0:(len(protein_length_nucleotides)-1)]#same for the variables counting stuff 
		continue
	x = line_list[0][len(line_list[0])-3:] #if this doesnt end in a stop codon then filter it out 
	stop_codon = 0
	if x == 'TAA': stop_codon += 1 
	if x == 'TGA': stop_codon += 1 
	if x == 'TAG': stop_codon += 1 
	if stop_codon == 0: #if the gene does not end with a stop codon then remove it from the list and count it as filterred out 
		stop_codon += 1
		genes = genes[0:(len(genes)-1)]
		ADAT_Codon_counts = ADAT_Codon_counts[0:len(ADAT_Codon_counts)-1] #remove via slcing again
		protein_length_nucleotides = protein_length_nucleotides[0:(len(protein_length_nucleotides)-1)]
		filter_counter += 1 
		continue
	#filterring by multiples of three is not necessary

	######Block1######
	#this block counts ADAT-dependent codons
	counter = 0
	for i in range(int((len(line_list[0]))/3)):
		current_codon = line_list[0][i*3:(i*3)+3] 
		y = current_codon in codon_set
		if y == True: ADAT_Codon_counts[len(ADAT_Codon_counts)-1] += 1 
	for i in range(len(line_list[0])):
		protein_length_nucleotides[len(protein_length_nucleotides)-1] += 1 
	

	
	######Block2######
	#this block gets the percentage of codons in the protein that are adat dependent 
	total_codon_counter = 0
	tapslivr_codon_counter = 0
	for i in range(int((len(line_list[0]))/3)):
		total_codon_counter += 1
		current_codon = line_list[0][i*3:(i*3)+3] 
		y = current_codon in codon_set
		if y == True: tapslivr_codon_counter += 1 
	ADAT_Codon_freqs.append(tapslivr_codon_counter/total_codon_counter)

	######Block 2.1#########
	#Calculate adat-dependent freqs along gene bodies
	gene_body_dict[genes[-1]] = []
	counts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	perc_counts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	for i in range(int((len(line_list[0]))/3)):
		percentile = (i+1)/total_codon_counter
		if 0 < percentile <= 0.05: in_dex = 0
		if 0.05 < percentile <= 0.1: in_dex = 1
		if 0.1 < percentile <= 0.15: in_dex = 2
		if 0.15 < percentile <= 0.2: in_dex = 3
		if 0.2 < percentile <= 0.25: in_dex = 4
		if 0.25 < percentile <= 0.3: in_dex = 5
		if 0.3 < percentile <= 0.35: in_dex = 6
		if 0.35 < percentile <= 0.4: in_dex = 7
		if 0.4 < percentile <= 0.45: in_dex = 8
		if 0.45 < percentile <= 0.5: in_dex = 9
		if 0.5 < percentile <= 0.55: in_dex = 10
		if 0.55 < percentile <= 0.6: in_dex = 11
		if 0.6 < percentile <= 0.65: in_dex = 12
		if 0.65 < percentile <= 0.7: in_dex = 13
		if 0.7 < percentile <= 0.75: in_dex = 14
		if 0.75 < percentile <= 0.8: in_dex = 15
		if 0.8 < percentile <= 0.85: in_dex = 16
		if 0.85 < percentile <= 0.9: in_dex = 17
		if 0.9 < percentile <= 0.95: in_dex = 18
		if 0.95 < percentile <= 1: in_dex = 19

		current_codon = line_list[0][i*3:(i*3)+3] 
		y = current_codon in codon_set
		if y == True: counts[in_dex] += 1 
		perc_counts[in_dex] += 1 
	for i in range(len(counts)):
		gene_body_dict[genes[-1]].append((counts[i]/perc_counts[i]))

	
	








	######Block3######
	#get percentages of TAPSLIVR and normalize by background AA content 
	#test for significance via simple chi square 
	total_codon_counter = 0 #observed and expected 
	tapslivr_codon_counter = 0 #observed 
	tcount = 0
	acount = 0
	pcount = 0
	scount = 0
	lcount = 0
	icount = 0
	vcount = 0
	rcount = 0
	ccount = 0
	dcount = 0
	ecount = 0
	fcount = 0
	gcount = 0
	hcount = 0
	kcount = 0
	ncount = 0
	qcount = 0
	ycount = 0



	for i in range(int((len(line_list[0]))/3)):
		total_codon_counter += 1
		current_codon = line_list[0][i*3:(i*3)+3] 
		y = current_codon in codon_set
		if y == True: tapslivr_codon_counter += 1 
		detect_N = current_codon in translatedict
		if detect_N == False: continue 
		currentaa = translatedict[current_codon]
		if currentaa == ['T']: tcount += 1 
		if currentaa == ['A']: acount += 1 
		if currentaa == ['P']: pcount += 1 
		if currentaa == ['S']: scount += 1 
		if currentaa == ['L']: lcount += 1 
		if currentaa == ['I']: icount += 1 
		if currentaa == ['V']: vcount += 1 
		if currentaa == ['R']: rcount += 1 
		if currentaa == ['C']: ccount += 1 
		if currentaa == ['D']: dcount += 1 
		if currentaa == ['E']: ecount += 1 
		if currentaa == ['F']: fcount += 1 
		if currentaa == ['G']: gcount += 1 
		if currentaa == ['H']: hcount += 1 
		if currentaa == ['K']: kcount += 1 
		if currentaa == ['N']: ncount += 1 
		if currentaa == ['Q']: qcount += 1 
		if currentaa == ['Y']: ycount += 1 


	chisquareexp = np.array([0,0])
	chisquareobs = np.array([0,0])
	chisquareobs[0] = tapslivr_codon_counter
	chisquareobs[1] = total_codon_counter
	######Calculate the expected number of adat dependent codons##########################
	chisquareexp[0] = round((tcount*tfrac)+(acount*afrac)+(pcount*pfrac)+(scount*sfrac)+(lcount*lfrac)+(icount*ifrac)+(vcount*vfrac)+(rcount*rfrac)+(ccount*cfrac)+(dcount*dfrac)+(ecount*efrac)+(fcount*ffrac)+(gcount*gfrac)+(hcount*hfrac)+(kcount*kfrac)+(ncount*nfrac)+(qcount*qfrac)+(ycount*yfrac)) #calculate expected number based on TAPSIVR content and genomic background 
	######Calculate the expected number of adat dependent codons##########################
	chisquareexp[1] = total_codon_counter
	

	adat_exp.append(chisquareexp[0])
	adat_obs.append(tapslivr_codon_counter)
	protein_length_AA.append(total_codon_counter) 

	from scipy.stats import chisquare
	x = (chisquare(chisquareobs, f_exp=chisquareexp))
	if chisquareexp[0] > chisquareobs[0]: 
		ADAT_pvalues_t.append(1-x[1]) #transform p values such that they can be ordered by relative statistical enrichment 
	if chisquareexp[0] <= chisquareobs[0]: 
		ADAT_pvalues_t.append(-1+x[1])  
	ADAT_pvalues.append(x[1])  	
	bonferroni_corrected_chisquare.append(x[1]*zzz)
	
	if correction == 'bonferroni':
		if x[1]*zzz < 0.05 and chisquareexp[0] > chisquareobs[0]: 
			bonferroni_enriched_depleted.append('depleted')
	
		if x[1]*zzz < 0.05 and chisquareexp[0] <= chisquareobs[0]:
			bonferroni_enriched_depleted.append('enriched')
		if x[1]*zzz >= 0.05: bonferroni_enriched_depleted.append('nonsignificant')
		if chisquareexp[0] == 0:
			if chisquareobs[0] == 0:
				bonferroni_enriched_depleted.append('nonsignificant')





hexome.close()

if correction == 'benjamani-hochberg':
	#place P values in ascending order
	
	ordered_pvalues = sorted(ADAT_pvalues)
	original_order = []
	x = 1 
	for i in range(len(ordered_pvalues)):
		original_order.append(x)
		x+=1

	list1 = ADAT_pvalues
	list4 = ordered_pvalues
	list2 = original_order
	list3 = ADAT_pvalues_t

	zipped_lists = zip(list1, list2, list3, list4 )
	sorted_pairs = sorted(zipped_lists)

	tuples = zip(*sorted_pairs)
	list1, list2, list3, list4 = [ list(tuple) for tuple in  tuples]

	#print(list2)

	BH_criticalvals = []
	for i in range(len(list1)):
		BH_criticalvals.append(((i+1)/len(genes))*FDR)

	#print(len(BH_criticalvals))
	
	BH_enriched_depleted = []
	for i in range(len(BH_criticalvals)):
		if list1[i] < BH_criticalvals[i]:
			if list3[i] > 0: BH_enriched_depleted.append('depleted')
			
			if list3[i] < 0: BH_enriched_depleted.append('enriched')
		else: BH_enriched_depleted.append('nonsignificant')






	zipped_lists = zip(list2, list1,list3,BH_enriched_depleted, list4)
	sorted_pairs = sorted(zipped_lists)

	tuples = zip(*sorted_pairs)
	list2, list1, list3, BH_enriched_depleted, list4 = [ list(tuple) for tuple in  tuples]

	#print(list2)






if correction == 'bonferroni':
	g = open(sys.argv[3], 'w')
	sys.stdout = g
	print('Gene,adatcodoncount,adatcodonfreq,proteinlennucleotides,ADAT_transformed_ordered_pvalues,ADAT_pvalues_raw,ADAT_pvalues_corrected,chisquareadatenrichedstatus,expected_adat_codons,observed_adat_codons,protein_length_AA,feature')
	for i in range(len(genes)):
		print(genes[i]+','+str(ADAT_Codon_counts[i])+','+str(ADAT_Codon_freqs[i])+','+str(protein_length_nucleotides[i])+','+str(ADAT_pvalues_t[i])+','+str(ADAT_pvalues[i])+','+str(bonferroni_corrected_chisquare[i])+','+str(bonferroni_enriched_depleted[i])+','+str(adat_exp[i])+','+str(adat_obs[i])+','+str(protein_length_AA[i])+',gene')





if correction == 'benjamani-hochberg':
	g = open(sys.argv[3], 'w')
	sys.stdout = g
	print('Gene,adatcodoncount,adatcodonfreq,proteinlennucleotides,ADAT_transformed_ordered_pvalues,ADAT_pvalues_raw,chisquareadatenrichedstatus,expected_adat_codons,observed_adat_codons,protein_length_AA,feature')
	for i in range(len(genes)):
		print(genes[i]+','+str(ADAT_Codon_counts[i])+','+str(ADAT_Codon_freqs[i])+','+str(protein_length_nucleotides[i])+','+str(ADAT_pvalues_t[i])+','+str(ADAT_pvalues[i])+','+str(BH_enriched_depleted[i])+','+str(adat_exp[i])+','+str(adat_obs[i])+','+str(protein_length_AA[i])+',gene')












sys.stdout = original_stdout 

if correction == 'benjamani-hochberg':
	gene_enriched_dict = {}
	for i in range(len(genes)):
		gene_enriched_dict[genes[i]] = BH_enriched_depleted[i]

if correction == 'bonferroni':
	gene_enriched_dict = {}
	for i in range(len(genes)):
		gene_enriched_dict[genes[i]] = bonferroni_enriched_depleted[i]


g = open(sys.argv[4], 'w')
sys.stdout = g

print('Gene,perc_len,perc_adat,chisquareadatenrichedstatus')
for k,v in gene_body_dict.items():
	for i in range(len(v)):
		print(k+','+str((i+1)*5)+','+str(v[i])+','+gene_enriched_dict[k])





