# ADATscan
Scripts, results, and test files for running ADATscan

#File 1:
combined_ADATscan_results.csv
contains the combined results for all four species with a column labeling the entry with the species (used for plotting)

#File 2: 
ADATscan_dmel.py
ADATscan script tailored to the fasta file provided on flybase 
Usage: python ADATscan_dmel.py <background model filename> <results filename>
  
#File 3
dmel_ADATscan_results.csv
results obtained in this study when running ADATscan on the drosophila melanogaster exome 
  
#File 4
dmel_background_model.csv
background model for ADAT-dependent codon usage in drosophila melanogaster 
  
#File 5
dmel_cds.fasta
dmel exome file retrieved from flybase sept 16,2022 http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/ (newline characters removed)

#File 6
dmel_panther_analysis_fishers_exact_bonferroni.txt
PANTHER results for drosophila melanogaster derived from http://geneontology.org/

#File 7
ADATscan_CCDS.py
ADATscan script tailored to the fasta files provided on CCDS (human and mouse) 
Usage: python ADATscan_CCDS.py <background model filename> <results filename>

#File 8
gProfiler_hsapiens_9-16-2022_11-25-39 AM.csv
gene name file for human data derived from https://biit.cs.ut.ee/gprofiler/convert. Added columns for unique CCDS names and unique gene names 

#File 9
human_ADATscan_results.csv
results obtained in this study when running ADATscan on the human exome 

#File 10 
human_background_model.csv
background model for ADAT-dependent codon usage in humans

#File 11
human_cds.fasta
human exome file retrieved from CCDS sept 16,2022 https://ftp.ncbi.nlm.nih.gov/pub/CCDS/ (newline characters removed)

#File 12
human_panther_analysis_fishers_exact_bonferroni.txt
PANTHER results for humans derived from http://geneontology.org/

#File 13
gProfiler_mmusculus_9-16-2022_11-48-46 AM.csv
gene name file for mouse data derived from https://biit.cs.ut.ee/gprofiler/convert. Added columns for unique CCDS names and unique gene names 

#File 14
mouse_ADATscan_results.csv
results obtained in this study when running ADATscan on the mouse exome 

#File 15
mouse_background_model.csv
background model for ADAT-dependent codon usage in mice

#File 16
mouse_cds.fasta
mouse exome file retrieved from CCDS sept 16,2022 https://ftp.ncbi.nlm.nih.gov/pub/CCDS/ (newline characters removed)

#File 17
mouse_panther_analysis_fishers_exact_bonferroni.txt
PANTHER results for mice derived from http://geneontology.org/

#File 18
ADATscan_yeast.py
ADATscan script tailored to the fasta files provided on SGD (Scer) 
Usage: python ADATscan_yeast.py <background model filename> <results filename>

#File 19
yeast_ADATscan_results.csv
results obtained in this study when running ADATscan on the Scer exome 

#File 20
yeast_background_model.csv
background model for ADAT-dependent codon usage in Scer

#File 21
yeast_cds.fasta
Scer exome file retrieved from SGD sept 16,2022  (newline characters removed)
http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/

#File 22 
scer_panther_analysis_fishers_exact_bonferroni.txt
PANTHER results for Scer derived from http://geneontology.org/

#File 23 
ADATscan_generalized.py
Usage: python ADATscan_generalized.py <background model filename> <results filename> <input codons 1 at a time then type "n">

