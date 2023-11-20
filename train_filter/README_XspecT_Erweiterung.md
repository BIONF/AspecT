# XspecT-Erweiterung

Expands XspecT, so new filter for a genus can automatically be trained. It's main 
script is XspecT_trainer.py. The rest of the scripts are inside the python module 
train_filter. 

## Training new filter

XspecT_trainer.py uses command line arguments. The examples for using XspecT_trainer.py
are using Salmonella since this genus only has two defined species in the NCBI 
databases.

### Jellyfish

The program jellyfish is used to count distinct k-meres in the assemblies. For XspecT_
trainer.py to work jellyfish needs to be installed. It can be installed using bioconda:

`
conda install -c bioconda jellyfish
`

### Training examples

New filters with assemblies from NCBI RefSeq can be trained with the following line. The 
python libraries from [requirements.txt](..%2Frequirements.txt) need to be installed.

`
python XspecT_trainer.py Salmonella 1
`

Training filters with custom data can be done using the following line.

`
python XspecT_trainer.py Salmonella 2 -bf /path/to/concate_assemblies -svm 
/path/to/assemblies
`

All command line arguments are explained using the following line.

`
python XspecT_trainer.py -h
`

# Explanation of the scripts

## backup_filter.py

Creates a backup of all files needed for the species assignment by XspecT for a specific
genus. The backup will be done, if new filters will be created for a genus which 
already has trained filters.

## create_svm.py

Downloads the needed assemblies and trains a support-vector-machine for the genus.

## extract_and_concatenate.py

Unzips the downloaded assemblies. Concatenates assemblies per species that will be used 
to train the bloomfilters.

## get_paths.py

Functions that get specific paths.

## html_scrap.py

Updates a list of all NCBI RefSeq assembly accessions that have a taxonomy check result
of OK. The taxonomy check from NCBI RefSeq uses the ANI (average-nucleotide-
identity) to compute a result.

## interface_XspecT.py

Mostly functions that train new bloomfilters automatically. The functions were 
originally writen for XspecT in a non-automatic way and were updated.

## k_mer_count.py

Uses jellyfish to count distinct k-meres in every concatenated assembly. The highest
count will be used to compute the size of the bloomfilters.

## ncbi_api

A module which makes requests to the NCBI Datasets API.

### download_assemblies.py

The specific function that downloads assemblies from NCBI RefSeq using NCBI 
datasets.

### ncbi_assembly_metadata.py

Takes a dictionary with species and their taxon ID and asks NCBI for assemblies of
the species. Saves the collected accessions of the found and selected assemblies.

### ncbi_children_tree.py

Takes the name or ID of a genus and gives a list with all its species.

### ncbi_taxon_metadata.py

Takes a list with taxon and collects metadata like their scientific name and rank.


















