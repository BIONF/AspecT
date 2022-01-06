import BF_v2
from multiprocessing import Process, Pipe
import pickle
import glob
import os
from collections import Counter


def get_added_genomes():
    """ Reads in pickled list, returns none if no new genomes have been added"""
    with open(r'filter/FilterClonetypes.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)

    # IC1 to IC8 are not deletable. That means if the IC-list is not longer than 8, than
    # there are no new IC's
    if len(clonetypes) == 8:
        added = [None]
    else:
        # gives all added genomes after IC8
        added = clonetypes[8:]

    return added


def read_search(IC_lookup, reads, quick, pipe=None):
    with open(r'filter/FilterClonetypes.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    # initialising filter with database parameters
    BF = BF_v2.AbaumanniiBloomfilter(123000000)
    BF.set_arraysize(123000000)
    BF.set_hashes(7)
    BF.set_k(20)

    # Array Size 22.000.000
    paths = [r'filter/IC1.txt',
             r'filter/IC2.txt',
             r'filter/IC3.txt',
             r'filter/IC4.txt',
             r'filter/IC5.txt',
             r'filter/IC6.txt',
             r'filter/IC7.txt',
             r'filter/IC8.txt']

    if IC_lookup[8]:
        # added Genomes
        # IC1 to IC8
        # Selecting wanted slices
        for i in [7, 6, 5, 4, 3, 2, 1, 0]:
            if IC_lookup[i]:
                pass
            else:
                del clonetypes[i]
                del paths[i]

        # getting all added files
        temp = glob.glob('filter/added/*.txt')
        added = []
        if len(temp) == 0:
            pass
        else:
            # these for-loops are needed for sorting the paths
            # so they match with the pickle-list order
            for i in range(len(clonetypes)):
                for j in range(len(temp)):
                    if clonetypes[i] in temp[j]:
                        added.append(temp[j])

            paths.extend(added)

        BF.read_clonetypes(paths, clonetypes)

    else:
        # Only IC1 to IC8
        # Selecting wanted slices
        clonetypes = clonetypes[:8]
        for i in [7, 6, 5, 4, 3, 2, 1, 0]:
            if IC_lookup[i]:
                pass
            else:
                del clonetypes[i]
                del paths[i]

    BF.read_clonetypes(paths, clonetypes)
    BF.lookup_txt(reads, quick)
    score = BF.get_score()
    hits = BF.get_hits_per_filter()
    names = BF.get_names()
    BF.cleanup()
    del BF

    if pipe is not None:
        pipe.send([score, names, hits])
        pipe.close()
    else:

        return score, names, hits


def pre_processing():
    "Preprocesses the Bloomfilter-Matrix when the program is launched"
    with open(r'filter/FilterSpecies.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    # initialising filter with database parameters
    # kmer20 = 115000000
    # kmer31 = 122000000
    BF = BF_v2.AbaumanniiBloomfilter(115000000)
    BF.set_arraysize(115000000)
    BF.set_hashes(7)
    BF.set_k(20)
    paths = sorted(os.listdir(r"filter/species/"))
    for i in range(len(paths)):
        paths[i] = r"filter/species/" + paths[i]
    BF.read_clonetypes(paths, clonetypes)
    return BF


def read_search_spec(reads, quick, BF):
    "Searches sequence-data in Bloomfilter and gets kmer-hits"
    del reads[0]
    BF.lookup_txt(reads, quick)
    score = BF.get_score()
    score_sum = 0
    for i in score:
        score_sum += i
    hits = BF.get_hits_per_filter()
    names = BF.get_names()
    return score, names, hits


def single_oxa(reads, ext, pipe=None):
    """Uses the Bloomfilter module to lookup the OXA-genes"""
    # getting filters
    paths = os.listdir(r"filter/OXAs/")
    oxas = []
    for i in paths:
        oxas.append(i[:-4])

    for i in range(len(paths)):
        paths[i] = r"filter/OXAs/" + paths[i]


    # initialising filter with database parameters
    BF = BF_v2.AbaumanniiBloomfilter(80000)
    BF.set_arraysize(80000)
    BF.set_clonetypes(len(paths))
    BF.set_hashes(7)
    BF.set_k(20)
    # User Options

    # reading single OXA filters
    BF.read_clonetypes(paths, oxas)

    # starting Bloomfilter process, depends on filetype
    BF.lookup_oxa(reads, ext)

    score = BF.get_oxa_score()
    BF.cleanup()
    del BF

    if pipe is not None:
        pipe.send([score, oxas])
        pipe.close()
    else:
        return score, oxas


def oxa_and_IC_multiprocessing(IC_lookup, reads, ext, quick):
    """ Uses Multiprocessing to lookup OXA genes and Clonetypes at the same time """
    # Sources:
    # https://docs.python.org/3/library/multiprocessing.html#sharing-state-between-processes
    # https://stackoverflow.com/questions/7207309/python-how-can-i-run-python-functions-in-parallel
    # using pipes to Transfer data between functions
    parent_ic, child_ic = Pipe()
    parent_oxa, child_oxa = Pipe()

    if ext == 'fq' or ext == 'fastq':
        reads_ct = reads[:2000]
    else:
        reads_ct = reads
    p1 = Process(target=read_search, args=(IC_lookup, reads_ct, quick, child_ic))
    p1.start()
    p2 = Process(target=single_oxa, args=(reads, ext, child_oxa))
    p2.start()
    p1.join()
    p2.join()

    # getting results back from pipes
    results_ic = parent_ic.recv()  # has scores and names
    results_oxa = parent_oxa.recv()  # has scores and names

    return results_ic[0], results_ic[1], results_ic[2], results_oxa[0], results_oxa[1]
