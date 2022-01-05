These folders contain all the test-data.

=== Folder Species-Assignment ===

This Folder contains the results (Score-Vectors) from the taxonomic Assignment of over 3000 Genomes using different SVM-Parameters and Bloom Filters.
File names explained:
T1 = Training-Set 1 with up to 4 Genomes/BF
T2 = Training-Set 2 with 1 Genome/BF
poly = polynomial kernel-function
rbf = radial basis-function
1.0 and 1.5 refer to the Regulatory Parameter C of the SVM

Further files with the name Test-Set_X contain the Score-Vectors of Tests with different amounts of used k-mers
x = only every xth k-mer of the input-data was used

The file Training_data_spec.csv contains the training-data (support vectors) of the SVM


=== Folder Strain typing ===

The Folder Strain-Typing contains the Score-vectors of both reference-datasets. Assembled Genomes and the
Core-Genome of A. baumannii have been tested.
Each dataset has been tested with the full-run (all k-mers of the Input-Sequence have been tested)
and the quick-run (only every 10th k-mer has been tested). The used refernce is emntioned in the file-name.
The Format is the following:

NCBI Accession (GCF_...), Score IC1, Score IC2, Score IC3, ..., Score IC8

The files with 'Reads as Input' are the results of tested non-overlapping read-samples.
The Format is the following:

Filename, number of tested reads, Score IC1, Score IC2, Score IC3, ..., Score IC8

If the filename contains '...15% of the nucleotides changed', then ~15% of the nucleotides
in a read have been changed in order to test sequencing-methods with high error rates.


=== Folder OXA-Search ===

The file 'OXA search in genomes' contains the results of the OXA-Search in genomes, the Format is the following:

NCBI Accession (GCF_...), Score OXA-23, Score OXA-24, Score OXA-51, Score OXA-58


The file 'OXA genes' contains the results of the OXA-Search used on other OXA-genes from the four families:

correct OXA-Familie of the gene, OXA-Gene name, NCBI Accession, Score OXA-23, Score OXA-24, Score OXA-51, Score OXA-58


The files 'OXA-search in assemblys' and 'OXA-search in reads' contain the results of OXA-Search in Sequence-reds and theiren 
assemblys. 250.000 reads per .fq-file have been used. The format is:

filename, Score OXA-23, Score OXA-24, Score OXA-51, Score OXA-58

