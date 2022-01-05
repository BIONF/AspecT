import BF_v2
import os
import csv
import pickle
import sys
from Bio import SeqIO, SeqRecord, Seq
import search_filter


def write_file(name):
    """Adds new species name to the species-list"""
    itemlist = pickle.load(open(r'filter/FilterSpecies.txt', 'rb'))
    itemlist.append(name)
    itemlist.sort(key = lambda x: x.split(",")[-1][:2])
    with open(r'filter/FilterSpecies.txt', 'wb') as fp:
        pickle.dump(itemlist, fp)


def train_Core():
    """trains (concatenated-)genomes into BF and saves them"""
    files = os.listdir(r'filter\new_species')
    for i in range(len(files) -1, -1, -1):
        if 'fna' in files[i] or 'fasta' in files[i]:
            continue
        else:
            del files[i]
    for i in range(len(files)):
        #set BF-parameters
        BF = BF_v2.AbaumanniiBloomfilter(115000000)
        BF.set_arraysize(115000000)
        BF.set_clonetypes(1)
        BF.set_hashes(7)
        BF.set_k(20)
        path = r'filter/new_species/' + files[i]
        name = files[i].split('.')[-2] + '.txt'
        result = r'filter/species/' + name
        BF.train_sequence(path, 0)
        BF.save_clonetypes(result)
        print("Added A. " + files[i].split('.')[-2] + " to AspecT")
        BF.cleanup()



def main():
    if len(sys.argv) != 1:
        print("Training Bloom Filters...")
        train_Core()
        for i in range(len(sys.argv)):
            if i == 0:
                continue
            # saving new species names
            write_file(sys.argv[i])
        print("Calculating new Training-Data...")
        BF = search_filter.pre_processing()
        # generate new training-data
        BF.helper()
        print("Finished!")
    else:
        print("Error: No Species-Name entered!")

if __name__ == '__main__':
    main()
