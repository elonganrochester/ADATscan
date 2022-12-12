#!/usr/bin/env python3
"""
python ADATscan_CCDS.py <fasta> <null model output> <Enrichment and Depletion output> 

<fasta> should lack newline characters within the genes and should only contain coding seq (no introns)

<output> is a file with 

Note: null model of codon usage is generated using the fasta as input

Note: if isoforms derived from the same coding sequence are not denoted with <gene>.1 <gene>.2 then modification is needed  

Note: identical identifiers are skipped and only the first instance is kept 
"""


#Import modules needed 
import sys
import numpy as np
import scipy
from scipy import stats

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

#ADAT dependent codons (no cognate trna in humans, no decoding without inosine)
codon_set = ['ACC','GCC','CCC','TCC','CTC','ATC','GTC','CGC']
##################STEP 1 GENERATE DATA STRUCTURES############################## 





##################STEP 2 GENERATE NULL MODEL OF CODON BIAS IN SAMPLE############################## 
#loop through and get genome wide codon bias for the background null model 
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
counter = 0 #initialize counter to count total lines in file 
#loop through exome file and extract codon usage 


start_count = 0
stop_count = 0
skipped = False 
for line in hexome:
	counter += 1 
	if skipped == True: #keeps from filtering out isoforms and entries that do not begin with start codon or end with stop codon twice, see below
		skipped = False
		continue	
	line_list = line.strip().split('|') #parse the line by the character following the gene name 
	duplicate = line_list[0][1:] in genes #skip duplicate identifiers in the file
	if duplicate == True: 
		skipped = True
		continue



	if line_list[0][0] == '>': # if the line is a header line record the gene name 
		if line_list[0][1:][-2:] != '.1':#filter out isoforms so same codons arent counted twice
			skipped = True 
			continue 
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
		currentaa = translatedict[current_codon] #translate the codon using dictionary made above 
		#count the TAPSLIVR aas and how many times the ADAT requiring codons are used 
		
		if currentaa == ['T']: tcount += 1 
		if current_codon == 'ACC' : tadatreqcount +=1 
		if currentaa == ['A']: acount += 1 
		if current_codon == 'GCC' : aadatreqcount +=1 
		if currentaa == ['P']: pcount += 1 
		if current_codon == 'CCC' : padatreqcount +=1 
		if currentaa == ['S']: scount += 1 
		if current_codon == 'TCC' : sadatreqcount +=1 
		if currentaa == ['L']: lcount += 1 
		if current_codon == 'CTC' : ladatreqcount +=1 
		if currentaa == ['I']: icount += 1 
		if current_codon == 'ATC' : iadatreqcount +=1 
		if currentaa == ['V']: vcount += 1 
		if current_codon == 'GTC' : vadatreqcount +=1 
		if currentaa == ['R']: rcount += 1 
		if current_codon == 'CGC' : radatreqcount +=1 
hexome.close()#close exome file 

#calculate genomic background of codon usage 
tfrac = tadatreqcount/tcount
afrac = aadatreqcount/acount
pfrac = padatreqcount/pcount
sfrac = sadatreqcount/scount
lfrac = ladatreqcount/lcount
ifrac = iadatreqcount/icount
vfrac = vadatreqcount/vcount
rfrac = radatreqcount/rcount

#Print general results of this upstream loop to the null model file 
original_stdout = sys.stdout 
f = open(sys.argv[2], 'w') #initialize strains, category, population, junction, homozygous vs heterozygous file (used to make donor file and hap file) 
sys.stdout = f
print("Genes lacking start codon count:", start_count,"Genes lacking stop codon count:", stop_count)
print("Genes retained from file:", len(genes))
print("Total genes in file:",int(counter/2))
print('ADAT fraction threonine:'+str(tfrac))
print('ADAT fraction alanine:'+str(afrac))
print('ADAT fraction proline:'+str(pfrac)) 
print('ADAT fraction serine:'+str(sfrac ))
print('ADAT fraction luecine:'+str(lfrac ))
print('ADAT fraction isoluecine:'+str(ifrac ))
print('ADAT fraction valine:'+str(vfrac ))
print('ADAT fraction arginine:'+str(rfrac)) 
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

protein_length_nucleotides = [] #double check to ensure the number of nucleotides is right for each gene
adat_exp = []
adat_obs = []
protein_length_AA = [] 
bonferroni_corrected_chisquare = []
bonferroni_enriched_depleted = []
skipped = False

start_count = 0
stop_count = 0
gene_body_dict = {}

for line in hexome: #for line in fasta file 
	if skipped == True: #for filtered out isoforms see below 
		skipped = False
		continue
	line_list = line.strip().split('|') #parse the line by the character following the gene name 
	if line_list[0][0] == '>': # if the line is a header line record the gene name 
		if line_list[0][1:][-2:] != '.1':
			skipped = True #skip the seq line after isoforms too above 
			continue 
		
		duplicate = line_list[0][1:] in genes
		if duplicate == True: 
			skipped = True
			continue

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
	for i in range(int((len(line_list[0]))/3)):
		total_codon_counter += 1
		current_codon = line_list[0][i*3:(i*3)+3] 
		y = current_codon in codon_set
		if y == True: tapslivr_codon_counter += 1 
		currentaa = translatedict[current_codon]
		if currentaa == ['T']: tcount += 1 
		if currentaa == ['A']: acount += 1 
		if currentaa == ['P']: pcount += 1 
		if currentaa == ['S']: scount += 1 
		if currentaa == ['L']: lcount += 1 
		if currentaa == ['I']: icount += 1 
		if currentaa == ['V']: vcount += 1 
		if currentaa == ['R']: rcount += 1 
	chisquareexp = np.array([0,0])
	chisquareobs = np.array([0,0])
	chisquareobs[0] = tapslivr_codon_counter
	chisquareobs[1] = total_codon_counter
	######Calculate the expected number of adat dependent codons##########################
	chisquareexp[0] = round((tcount*tfrac)+(acount*afrac)+(pcount*pfrac)+(scount*sfrac)+(lcount*lfrac)+(icount*ifrac)+(vcount*vfrac)+(rcount*rfrac)) #calculate expected number based on TAPSLIVR content and genomic background 
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





