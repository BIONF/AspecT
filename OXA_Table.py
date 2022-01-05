import json
from Bio import SeqIO
import os


class OXATable:

    def __init__(self):
        self.kmere = {}
        self.total = 0
        self.found = 0

    def create_table(self, directory):
        """ Reads in fasta file, creates Dictionary with k-mer counter"""
        # taking all fasta files for one big table
        oxas = {}
        files = os.listdir(directory)
        for i in range(len(files)):
            kmere = {}
            file = directory + files[i]
            for sequence in SeqIO.parse(file, "fasta"):
                for j in range(0, len(sequence.seq) - 20 + 1):
                    kmer = str(sequence.seq[j: j + 20])
                    count = kmere.get(kmer, 0)
                    kmere[kmer] = count + 1
            oxas[files[i][:-6]] = kmere

    def lookup(self, gene, kmer):
        """ Tests if kmer in dictionary, if so: reduces the counter"""

        # Only returns True if kmer has been found and there was one left
        if kmer in self.kmere[gene]:
            if self.kmere[gene][kmer] > 0:
                self.kmere[gene][kmer] -= 1
                return True
            else:
                return False
        else:
            return False

    def save_dic(self, path):
        """ writes dictionary to file using json"""
        json.dump(self.kmere, open(path, 'w'))

    def read_dic(self, path):
        """ Reads dictionary from file using json"""
        self.kmere = json.load(open(path))

    def cleanup(self):
        self.kmere = {}

    def get_counter(self, path=r'filter/OXAs_dict/counter.txt'):
        counter = json.load(open(path))
        return counter

