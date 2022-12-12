# ADATscan
This repository contains the code for ADATscan and results derived from ADATscan. 
ADATscan is a python script that takes in a fasta file of orfs (no introns, no newline characters in sequences)
ADATscan requires the user to input the multiple testing correction method they desire. Bonferroni or the Banjamini-Hochberg have been implemented. If the Banjamini-Hochberg correction is desired then the user must also input the desired FDR.  
ADATscan outputs three files: 
1) The background model of ADAT-dependent codon usage in the exome provided
2) The results of the statistical tests and summary stats for each gene 
3) A file tyhat allows for plotting ADAT-dependent codon freqeucnies along gene bodies among enriched, depleted, and nonsignificant genes.
 
There are five versions of the script in this repository
1) Generalized version
2) tailored for CCDS
3) tailored for wormbase 
4) tailored for flybase 
5) tailored for sgd 

Notation in these databases are idiosyncratic, hence serarate versions. When using the generalized version, it is VERY important that redundent orfs (for example isoforms) are excluded. The background model will be calculated incorrectly if isoforms are retained. 

The example folder contains a word document that walks through downloading and running ADATscan.
The scripts folder contains a 
The results folder contains data derived from ADATscan for human, mouse, nematode, fruit fly, and yeast exomes. Data derived from analyses using the Bonferroni correction and the Benjamini-Hochberg procedure (FDR = 0.01) are prvided in separate directories. 

When using the software provided as it pertains to ADAT-dependent codons, it is critical that the subset of focal codons is correct. Check appropriate databases for tRNA genes that are present in the species of interest. 

License: See LICENSE.txt (MIT) 


