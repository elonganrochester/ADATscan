# ADATscan
Scripts and sample results for running ADATscan

#File 1:
combined_ADATscan_results.csv
contains the combined results for all four species with a column labeling the entry with the species (used for plotting)

#File 2: 
ADATscan_dmel.py
ADATscan script tailored to the fasta file provided on flybase 
Usage: python ADATscan_dmel.py <input fasta exome file (no newlines other than those seperating entry names from entries!)> <background model filename> <results filename>
  
#File 3
dmel_ADATscan_results.csv
results obtained in this study when running ADATscan on the drosophila melanogaster exome 
  
#File 4
dmel_background_model.csv
background model for ADAT-dependent codon usage in drosophila melanogaster 
  
#File 5
dmel_panther_analysis_fishers_exact_bonferroni.txt
PANTHER results for drosophila melanogaster derived from http://geneontology.org/

#File 6
ADATscan_CCDS.py
ADATscan script tailored to the fasta files provided on CCDS (human and mouse) 
Usage: python ADATscan_CCDS.py <input fasta exome file (no newlines other than those seperating entry names from entries!)> <background model filename> <results filename>

#File 7
gProfiler_hsapiens_9-16-2022_11-25-39 AM.csv
gene name file for human data derived from https://biit.cs.ut.ee/gprofiler/convert. Added columns for unique CCDS names and unique gene names 

#File 8
human_ADATscan_results.csv
results obtained in this study when running ADATscan on the human exome 

#File 9
human_background_model.csv
background model for ADAT-dependent codon usage in humans

#File 10
human_panther_analysis_fishers_exact_bonferroni.txt
PANTHER results for humans derived from http://geneontology.org/

#File 11
gProfiler_mmusculus_9-16-2022_11-48-46 AM.csv
gene name file for mouse data derived from https://biit.cs.ut.ee/gprofiler/convert. Added columns for unique CCDS names and unique gene names 

#File 12
mouse_ADATscan_results.csv
results obtained in this study when running ADATscan on the mouse exome 

#File 13
mouse_background_model.csv
background model for ADAT-dependent codon usage in mice

#File 14
mouse_panther_analysis_fishers_exact_bonferroni.txt
PANTHER results for mice derived from http://geneontology.org/

#File 15
ADATscan_yeast.py
ADATscan script tailored to the fasta files provided on SGD (Scer) 
Usage: python ADATscan_yeast.py <input fasta exome file (no newlines other than those seperating entry names from entries!)> <background model filename> <results filename>

#File 16
yeast_ADATscan_results.csv
results obtained in this study when running ADATscan on the Scer exome 

#File 17
yeast_background_model.csv
background model for ADAT-dependent codon usage in Scer

#File 18 
scer_panther_analysis_fishers_exact_bonferroni.txt
PANTHER results for Scer derived from http://geneontology.org/

#File 19 
ADATscan_generalized.py
Usage: python ADATscan_generalized.py <input fasta exome file (no newlines other than those seperating entry names from entries!)> <background model filename> <results filename> <input codons 1 at a time then type "n">
NOTE: ALTERNATIVE GENETIC CODES MUST BE SPECIFIED IN THE SCRIPT IF RUNNING ON AN EXOME OF A SPECIES WITH AN ALTERRED GENETIC CODE

#File 20 
ADATscan_celegans.py
ADATscan script tailored to the fasta files provided on Wormbase
Usage: python ADATscan_celegans.py <input fasta exome file (no newlines other than those seperating entry names from entries!)> <background model filename> <results filename>

#File 21
c_elegans_adatscan_background_model.csv
background model for ADAT-dependent codon usage in C elegans

#File 22 
c_elegans_adatscan_results.csv
results obtained in this study when running ADATscan on the C elegans exome 

