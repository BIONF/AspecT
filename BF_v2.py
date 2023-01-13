try:
    # try with a fast c-implementation ...
    import mmh3 as mmh3
except ImportError:
    # ... otherwise fallback to this code!
    import pymmh3 as mmh3
from bitarray import bitarray
from Bio import SeqIO
from copy import deepcopy
from OXA_Table import OXATable
from Bio.Seq import Seq
import os
import csv
import time



class AbaumanniiBloomfilter:
    """ Bloomfilter that can read FASTA and FASTQ files to assign the given file to a reference-genome"""
    # Implementation of the Bloomfilter Project for Acinetobacter baumannii
    # Used an customized also for the Bloomfilter Project for Acinetobacter Species Assignment
    # Variables from the Strain-Typing were used if possible for the Species-Assignment to not over-complicate the Code
    # Code partly from https://github.com/Phelimb/BIGSI

    clonetypes = 1  # Number of IC's/Species
    hits_per_filter = [0] * clonetypes  # Hit counter per IC/per Species
    array_size = 22000000  # Standard arraysize per IC is 22mio for Core-genome
    hashes = 7  # Number of used Hash-functions
    k = 20  # length of the k-meres
    names = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'IC6', 'IC7', 'IC8']  # names of the IC's
    number_of_kmeres = 0  # counter of k-meres, will be used to calculate score
    reads = 1000  # standard read number

    def __init__(self, arraysize):
        """ creates empty matrix"""
        pass
        self.matrix = bitarray(arraysize)
        self.matrix.setall(False)
        self.array_size = arraysize
        self.kmeres = []
        self.hits_per_filter_kmere = []
        self.coverage = []
        self.hit = False
    # Setter


    def set_arraysize(self, new):
        """ changes Arraysize to new input-value, does not recreate matrix"""
        self.array_size = new


    def set_clonetypes(self, new):
        """ changes number of Clonetypes"""
        self.clonetypes = new
        self.hits_per_filter = [0] * self.clonetypes


    def set_hashes(self, new):
        """Changes number of used hash-functions"""
        self.hashes = new


    def set_k(self, new):
        """ Changes length of k-meres"""
        self.k = new


    def set_names(self, new):
        """ Changes Names of Filters, Input must be a List of names"""
        self.names = new


    def reset_counter(self):
        """resets counter"""
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes


    def set_reads(self, new):
        """ Changes number of reads to new value"""
        self.reads = new

    # Getter

    def get_score(self):
        """calculates score for all clonetypes
            Score is #hits / #kmeres"""

        score = []

        # calculates float for each value in [hits per filter]
        for i in range(self.clonetypes):
            if self.hits_per_filter[i] == 0:
                score.append(0.0)
            else:
                score.append(round(float(self.hits_per_filter[i]) / float(self.number_of_kmeres), 2))

        return score


    def get_norm(self):
        """ Divides each vector entry by sum of vector entrys"""
        s = sum(self.hits_per_filter)
        score = []

        # calculates float for each value in [hits per filter]
        for i in range(self.clonetypes):
            if self.hits_per_filter[i] == 0 or s == 0:
                score.append(0.0)
            else:
                score.append(round(float(self.hits_per_filter[i]) / s, 2))

        return score


    def get_reads(self):
        """ gets number of reads """
        return self.reads


    def get_hits_per_filter(self):
        """gets Hits per Filter"""
        return self.hits_per_filter


    def get_kmeres_per_sequence(self):
        """gets K-mer counter"""
        # returns number of k-meres per file
        return self.number_of_kmeres


    def get_names(self):
        """ gets names of filters"""
        return self.names


    def get_coverage(self):
        """ gets coverage"""
        return self.coverage

    # File management

    def save_clonetypes(self, path):
        """saves matrix as a binary file to the input-path"""
        # saving filters of clonetypes

        # creating file and saving matrix with the bitarray modul
        with open(path, 'wb') as fh:
            # writing to file with bitarray command
            self.matrix.tofile(fh)


    def read_clonetypes(self, paths, names):
        """ reads slices from files and concats them to a matrix,
        paths is list of paths and names is a string list"""

        # Updating parameters
        self.clonetypes = len(paths)
        self.names = names
        self.matrix = bitarray(0)
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes
        i = 0

        # creating matrix from single filters
        for path in paths:

            temp = bitarray()

            with open(path, 'rb') as fh:
                temp.fromfile(fh)
            i += 1
            self.matrix.extend(temp)

    # Bloomfilter

    def hash(self, kmer):
        """Hashes given string and returns Positions for the Array"""

        # Empty list for Array positions
        positions = []
        # Creating hashes for needed number of hash functions
        for i in range(self.hashes):
            # mmh3 takes that string and a seed,
            # each hash function takes an individual seed
            # after that, the hash-value will me divided by the array size until
            # a position in the array is guaranteed
            positions.append(mmh3.hash(kmer, i) % self.array_size)

        return positions


    def lookup(self, kmer, limit=False):
        """checks if an element is in the filters, returns list with True/False,
           takes kmer input string and checks all clonetypes if the k-mer is inside that set of kmers"""

        # getting positions
        positions = self.hash(str(kmer))
        # control if element is in filter
        hits = [True] * self.clonetypes
        self.hit = False

        for i in range(self.clonetypes):
            row = i*self.array_size
            # all 7 Positions are hardcoded, the number of hashes is always(!) 7
            # if all positions  are True, then hits[i] will also stay True
            # (i*self.array_size) skips to the same position in the next filter
            hits[i] = (self.matrix[positions[0] + row] &
                       self.matrix[positions[1] + row] &
                       self.matrix[positions[2] + row] &
                       self.matrix[positions[3] + row] &
                       self.matrix[positions[4] + row] &
                       self.matrix[positions[5] + row] &
                       self.matrix[positions[6] + row])

            if hits[i]:
                self.hit = True
                if limit:
                    if self.table.lookup(self.names[i], kmer):
                        self.hits_per_filter[i] += 1
                else:
                    # Update hit counter
                    self.hits_per_filter[i] += 1


    def train(self, kmer, clonetype):
        """ trains specific filter for a k-mer, input is that kmer and the desired Filter"""

        # getting hash Values
        positions = self.hash(kmer)
        # changing 0s to 1 in filter
        for i in range(len(positions)):
            # getting position of cell
            self.matrix[self.array_size * clonetype + positions[i]] = True


    def train_sequence(self, filepath, clonetype, quick=False):
        """trains whole sequence into filter, takes filepath to file and the desired filter as input"""
        # for each sequence (in multi-FASTA file)
        if quick:
            for sequence in SeqIO.parse(filepath, "fasta"):
                # for each k-mere
                for i in range(len(sequence.seq) - self.k):
                    # trains k-mere into filter
                    self.train(str(sequence.seq[i: i + self.k]), clonetype)
        else:
            for sequence in SeqIO.parse(filepath, "fasta"):
                # for each k-mere
                # for i in range(len(sequence.seq) - self.k + 1):
                for i in range(len(sequence.seq) - self.k + 1):
                    # trains k-mere into filter
                    self.train(str(sequence.seq[i: i + self.k]), clonetype)
                    # testing
                    self.train(str(sequence.seq[i: i + self.k].reverse_complement()), clonetype)



    def train_lines(self, lines, ct):
        """ Trains Extracted lines of fasta/fna file, given as list of strings"""
        for j in range(len(lines)):
            for i in range(len(lines[j]) - self.k + 1):
                # trains k-mere into filter
                self.train(str(lines[j][i: i + self.k]), ct)


    def lookup_sequence(self, path):
        """uses lookup function for whole sequence, takes path to file: file must be FASTA"""

        # Counter of k-meres
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes

        # for each sequence (in multi-FASTA file)

        # counting hits in Intervals
        self.readset = [0] * self.clonetypes

        for sequence in SeqIO.parse(path, "fasta"):

            # for each k-mere
            for i in range(0, len(sequence.seq) - self.k + 1):

                # lookup for all k-meres in filter
                self.lookup(str(sequence.seq[i: i + self.k]))
                self.number_of_kmeres += 1


    def lookup_txt(self, reads, ext, quick=False):
        """ Reading extracted fq-reads"""
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes

        if quick == 1:
            # Quick: Non-overlapping k-mers
            # AspecT-Quick-Mode every 500th kmer
            for single_read in reads:
                # r is rest, so all kmers have size k
                for j in range(0, len(single_read) - self.k, 500):
                    if "N" in single_read[j:j + self.k]:
                        continue
                    self.number_of_kmeres += 1
                    self.lookup(single_read[j: j + self.k])
        # AspecT Sequence-Reads every 10th kmer
        elif quick == 2:
            for single_read in range(0, len(reads)):
                hit_counter = 0
                for j in range(0, len(reads[single_read]) - self.k, 10):
                    if j == 5 and hit_counter == 0:
                        break
                    # updating counter
                    self.number_of_kmeres += 1
                    # lookup for kmer
                    temp = reads[single_read]
                    self.lookup(str(temp[j: j + self.k]))
                    if self.hit == True:
                        hit_counter += 1
        elif quick == 3:
            #ClAssT Quick-Mode every 10th kmer
            for single_read in reads:
                # r is rest, so all kmers have size k
                for j in range(0, len(single_read) - self.k, 10):
                    if "N" in single_read[j:j + self.k]:
                        continue
                    self.number_of_kmeres += 1
                    self.lookup(single_read[j: j + self.k])
        # metagenome mode
        elif quick == 4:
            counter = 0
            #print("Kmer Anzahl: ", len(reads))
            for kmer in reads:
                counter += 1
                # lookup for kmer
                hits_per_filter_copy = self.hits_per_filter[:]
                self.lookup(kmer)
                if ((sum(self.hits_per_filter) - sum(hits_per_filter_copy)) <= 5) and ((sum(self.hits_per_filter) - sum(hits_per_filter_copy)) != 0):
                    self.number_of_kmeres += 1
                elif (sum(self.hits_per_filter) - sum(hits_per_filter_copy)) > 5:
                    self.hits_per_filter = hits_per_filter_copy[:]
                if ext == "fasta" or ext == "fna" or ext == "fa":
                    if counter >= 50000:
                        break
        else:
            for single_read in reads:
                for j in range(len(single_read) - self.k + 1):
                    # updating counter
                    self.number_of_kmeres += 1
                    # lookup for kmer
                    self.lookup(str(single_read[j: j + self.k]))


    def lookup_unique(self, reads, quick=False):
        """ Reading extracted fq-reads, only unique kmers"""
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes
        unique = set()
        if quick:
            # Quick: Non-overlapping k-mers
            for single_read in reads:
                # r is rest, so all kmers have size k
                for j in range(0, len(single_read) - self.k, 100):
                    unique.add(single_read[j: j + self.k])
            for kmer in unique:
                self.number_of_kmeres += 1
                self.lookup(kmer)
        else:
            for single_read in reads:
                for j in range(len(single_read) - self.k + 1):
                    # updating counter
                    self.number_of_kmeres += 1
                    # lookup for kmer
                    self.lookup(str(single_read[j: j + self.k]))


    def helper(self):
        """Creates svm Traingings-Data from a set of genomes"""
        #https://pythonguides.com/python-write-a-list-to-csv/
        #https://stackoverflow.com/questions/21431052/sort-list-of-strings-by-a-part-of-the-string
        files = os.listdir(r'Training_data\genomes')
        # delete all non fna/fasta files from list
        for i in range(len(files) -1, -1, -1):
            if 'fna' in files[i] or 'fasta' in files[i]:
                continue
            else:
                del files[i]
        paths = files[:]
        scores = []
        files_split = []
        names = []
        #extracts the GCF-Number
        for i in range(len(files)):
            files_split.append(files[i].split("_"))
            try:
                files_split[i] = files_split[i][0] + "_" + files_split[i][1]
            except:
                files_split[i] = files[i].split('.')[-2]
        for i in range(len(files)):
            paths[i] = r'Training_data/genomes/' + paths[i]
        #extracts the names of all species from the file-name
        for i in range(len(files)):
            with open(paths[i]) as file:
                head = file.readline()
                head = head.split()
                if head[2] == "sp.":
                    names.append("none")
                    continue
                names.append(head[2])
        #performs a lookup in the BF and saves the scores in a list
        for i in range(len(files)):
            self.number_of_kmeres = 0
            self.hits_per_filter = [0] * self.clonetypes
            for sequence in SeqIO.parse(paths[i], "fasta"):
                for j in range(0, len(sequence.seq) - self.k, 500):
                    self.number_of_kmeres += 1
                    self.lookup(str(sequence.seq[j: j + self.k]))
            score = self.get_score()
            score = [str(x) for x in score]
            score = ",".join(score)
            scores.append(files_split[i] + "," + score + "," + names[i])
        #sorts the list by species name
        scores.sort(key = lambda x: x.split(",")[-1][:2])
        names = [x for x in names if x != "none"]
        names = list(dict.fromkeys(names))
        scores.insert(0, sorted(names))
        scores[0] = ["File"] + scores[0] + ["Label"]
        for i in range(1, len(scores)):
            scores[i] = [scores[i]]
        #writes the Traingings-Data to a csv-filter
        with open(r'Training_data/Training_data_spec.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(scores)


    def cleanup(self):
        """deletes matrix"""
        self.matrix = None


    def lookup_oxa(self, reads, ext):
        """ Looks for OXA Genes: Extension (ext) selects the fq-seach or fasta-search mode"""
        self.table = OXATable()
        self.table.read_dic(r'filter/OXAs_dict/oxa_dict.txt')
        if ext == 'fq':
            # fq mode
            for i in range(len(reads)):
                # going through all reads, discarding those who don't get any hits with 3 test k-meres

                # Building 3 test-kmeres: first, last, and middle
                k1 = reads[i][0:self.k]  # first k-mer
                k2 = reads[i][len(reads[i]) - self.k:]  # last k-mer
                mid = len(reads[i])//2
                k3 = reads[i][mid:mid+self.k]  # k-mer in middle

                # Taking sum of list as reference, if sum has not increased after testing those 3 kmeres,
                # then the read won't be tested further
                hit_sum = sum(self.hits_per_filter)
                copy = deepcopy(self.hits_per_filter)
                self.lookup(k1, True)
                self.lookup(k2, True)
                self.lookup(k3, True)

                # needs at least 2 of 3 hits to continue with read
                if (sum(self.hits_per_filter) - hit_sum) > 1:

                    for j in range(1, len(reads[i]) - 1 - self.k + 1):
                        # Skipping first, last and middle k-mer
                        if j != mid:
                            self.lookup(reads[i][j:j + self.k], True)
                            self.number_of_kmeres += 1

                else:
                    # resetting hit counter
                    self.hits_per_filter = copy

                # same, but with reverse complement
                reads[i] = Seq(reads[i])
                reads[i] = reads[i].reverse_complement()
                k1 = reads[i][0:self.k]  # first k-mer
                k2 = reads[i][len(reads[i]) - self.k:]  # last k-mer
                mid = len(reads[i]) // 2
                k3 = reads[i][mid:mid + self.k]  # k-mer in middle

                # Taking sum of list as reference, if sum has not increased after testing those 3 kmeres,
                # then the read won't be tested further
                hit_sum = sum(self.hits_per_filter)
                copy = deepcopy(self.hits_per_filter)
                self.lookup(k1, True)
                self.lookup(k2, True)
                self.lookup(k3, True)

                # needs at least 2 of 3 hits to continue with read
                if (sum(self.hits_per_filter) - hit_sum) > 1:

                    for j in range(1, len(reads[i]) - 1 - self.k + 1):
                        # Skipping first, last and middle k-mer
                        if j != mid:
                            self.lookup(reads[i][j:j + self.k], True)
                            self.number_of_kmeres += 1

                else:
                    # resetting hit counter
                    self.hits_per_filter = copy

        else:
            # fasta mode
            # Altes testen mit Genom, hits per filter ausgeben lassen
            #self.oxa_search_genomes(reads)
            #self.oxa_search_genomes_v2(reads)
            coordinates_forward = self.oxa_search_genomes_v3(reads)
            reads_reversed = []
            for r in range(len(reads)):
                # building reverse complement
                reads_reversed.append(Seq(reads[r]))
                reads_reversed[r] = reads_reversed[r].reverse_complement()
            # lookup reverse complement
            #self.oxa_search_genomes(reads)
            #self.oxa_search_genomes_v2(reads)
            coordinates_reversed = self.oxa_search_genomes_v3(reads_reversed)

        # cleanup
        reads = None
        self.table.cleanup()
        return coordinates_forward, coordinates_reversed


    def oxa_search_genomes(self, genome):
        for i in genome:
            for j in range(0, len(i), 20):
                hits = sum(self.hits_per_filter)
                kmer = i[j:j+self.k]
                self.lookup(kmer, True)

                if sum(self.hits_per_filter) > hits:
                    for n in range(j - 19, j + 20, 1):
                        if 0 <= j < len(i):
                            kmer = i[n:n + self.k]
                            self.lookup(kmer, True)
                else:
                    pass


    def oxa_search_genomes_v2(self, genome):
        for i in genome:
            j = 0
            success = False
            while(j < len(i)):
                hits = sum(self.hits_per_filter)
                kmer = i[j:j+self.k]
                self.lookup(kmer, True)
                if success == False:
                    if sum(self.hits_per_filter) > hits:
                        for n in range(j - 19, j + 20, 1):
                            if 0 <= j < len(i):
                                kmer = i[n:n + self.k]
                                self.lookup(kmer, True)
                        j += 40
                        success = True
                    else:
                        j += 20
                        success = False
                else:
                    if sum(self.hits_per_filter) > hits:
                        for n in range(j, j + 20, 1):
                            if 0 <= j < len(i):
                                kmer = i[n:n + self.k]
                                self.lookup(kmer, True)
                        j += 40
                        success = True
                    else:
                        j += 20
                        success = False


    def oxa_search_genomes_v3(self, genome):
        coordinates = []
        for i in genome:
            j = 0
            success = False
            while(j < len(i)):
                hits = sum(self.hits_per_filter)
                kmer = i[j:j+self.k]
                self.lookup(kmer, True)
                if success == False:
                    if sum(self.hits_per_filter) > hits:
                        counter = 0
                        coordinates.append([j])
                        # 1024 (longest oxa-gene) - 19
                        for n in range(j - 249, j + 1005, 1):
                            if 0 <= j < len(i):
                                hits_per_filter_copy = self.hits_per_filter[:]
                                kmer = i[n:n + self.k]
                                self.lookup(kmer, True)
                                if hits_per_filter_copy != self.hits_per_filter:
                                    counter += 1
                        if counter > 300:
                            coordinates[-1].append(j+1005)
                        else:
                            coordinates.pop()
                        j += 1005
                        success = True
                    else:
                        #j += 20
                        j += 250
                        success = False
                else:
                    if sum(self.hits_per_filter) > hits:
                        coordinates.append([j])
                        counter = 0
                        for n in range(j, j + 1005, 1):
                            if 0 <= j < len(i):
                                kmer = i[n:n + self.k]
                                hits_per_filter_copy = self.hits_per_filter[:]
                                self.lookup(kmer, True)
                                if hits_per_filter_copy != self.hits_per_filter:
                                    counter += 1
                        if counter > 300:
                            coordinates[-1].append(j+1005)
                        else:
                            coordinates.pop()
                        j += 1005
                        success = True
                    else:
                        j += 250
                        success = False
        #if len(coordinates) > 0:
            #print("Coordinates: ", coordinates)
        return coordinates


    def get_oxa_score(self):
        """ Returning hits per OXA/kmere in OXA-filter"""
        table = OXATable()
        counter = table.get_counter()
        score = []
        # calculates float for each value in [hits per filter]
        for i in range(self.clonetypes):
            if self.hits_per_filter[i] == 0:
                score.append(0.0)
            else:
                score.append(round(float(self.hits_per_filter[i]) / float(counter[self.names[i]]), 2))
        self.hits_per_filter = [0] * self.clonetypes
        return score
