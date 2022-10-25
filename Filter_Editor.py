import BF_v2
import pickle
import csv
import os
import json
from Bio.Seq import Seq


def add_filter(name, svm, lines):
    """ Adds new filter to API """
    # Adding name to file that holds all names of added filters
    # https://wiki.python.org/moin/UsingPickle
    added = pickle.load(open(r'filter/FilterClonetypes.txt', "rb"))
    added.append(name)
    pickle.dump(added, open(r'filter/FilterClonetypes.txt', "wb"))

    # Changing SVM
    with open(r'Training_data/Training_data_IC.csv', "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(svm)

    # Adding file
    BF = BF_v2.AbaumanniiBloomfilter(123000000)
    BF.set_arraysize(123000000)
    BF.set_hashes(7)
    BF.set_k(20)
    BF.set_clonetypes(1)
    # Splitting because of different sections in fasta files
    lines = ''.join(lines)
    lines = lines.split('>')
    BF.train_lines(lines, 0)
    BF.save_clonetypes(r'filter/added/' + name + '.txt')

    # Cleanup
    BF.cleanup()
    BF = None


def remove_filter(name):
    """ Removes filter from API"""
    # Removing file
    os.remove(r'filter/added/' + name + '.txt')
    # removing name from pickle list
    added = pickle.load(open(r'filter/FilterClonetypes.txt', "rb"))
    index = added.index(name)
    del added[index]
    pickle.dump(added, open(r'filter/FilterClonetypes.txt', "wb"))

    # removing rows and cols from csv
    r = csv.reader(open(r'Training_data/Training_data_IC.csv', 'r'))
    svm = list(r)
    header = svm[0]
    svm = svm[1:]

    # Getting Index of Cols that should be deleted
    index = header.index(name)

    for i in range(len(svm)-1, -1, -1):

        if svm[i][-1] == name:
            # Deleting rows
            del svm[i]
        else:
            # Deleting Col
            del svm[i][index]

    # Putting modified header back in
    del header[index]
    svm.insert(0, header)

    # Putting new Data in Reference-file
    with open(r'Training_data/Training_data_IC.csv', "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(svm)


def edit_svm(svm):
    """ Changing SVM Training-vectors """
    # svm.insert(0, ['File', 'IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'IC6', 'IC7', 'IC8', 'Label'])
    with open(r'Training_data/Training_data_IC.csv', "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(svm)


def add_oxa(name, lines):
    # Creating Filter

    BF = BF_v2.AbaumanniiBloomfilter(80000)
    BF.set_arraysize(80000)
    BF.set_hashes(7)
    BF.set_k(20)
    BF.set_clonetypes(1)
    # Splitting because of different sections in fasta files
    lines = ''.join(lines)
    lines = lines.split('>')
    for i in range(len(lines) - 1, -1, -1):
        if lines[i] == '':
            del lines[i]
    lines_reverse = []
    for r in range(len(lines)):
        # building reverse complement
        lines_reverse.append(Seq(lines[r]))
        lines_reverse[r] = str(lines_reverse[r].reverse_complement())
    lines = lines_reverse
    BF.train_lines(lines, 0)
    BF.save_clonetypes(r'filter/OXAs/' + name + '.txt')

    # Cleanup
    BF.cleanup()
    BF = None

    # Adding kmeres
    kmere = {}
    for i in range(len(lines)):
        for j in range(len(lines[i]) - 20 + 1):
            kmer = str(lines[i][j: j + 20])
            count = kmere.get(kmer, 0)
            kmere[kmer] = count + 1

    file = json.load(open(r'filter/OXAs_dict/oxa_dict.txt'))
    file[name] = kmere
    json.dump(file, open(r'filter/OXAs_dict/oxa_dict.txt', 'w'))
    file = None

    # Adding counter
    c = 0
    for i in range(len(lines)):
        if len(lines[i]) < 20:
            c += 1
        else:
            c += len(lines[i]) - 20 + 1

    counter = json.load(open(r'filter/OXAs_dict/counter.txt'))
    counter[name] = c
    json.dump(counter, open(r'filter/OXAs_dict/counter.txt', 'w'))

    lines = None


def remove_oxa(name):
    # Deleting Filter
    os.remove(r'filter/OXAs/' + name + '.txt')

    # Deleting k-mer counter
    counter = json.load(open(r'filter/OXAs_dict/counter.txt'))
    del counter[name]
    json.dump(counter, open(r'filter/OXAs_dict/counter.txt', 'w'))
    counter = None

    # Deleting kmeres
    file = json.load(open(r'filter/OXAs_dict/oxa_dict.txt'))
    del file[name]
    json.dump(file, open(r'filter/OXAs_dict/oxa_dict.txt', 'w'))
    file = None
