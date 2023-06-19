import BF_v2
from multiprocessing import Process, Pipe
import pickle
import glob
import os
from collections import Counter
import time


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
    BF.lookup_txt(reads, False, quick)
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


def pre_processing_ClAssT():
    "Preprocesses the Bloomfilter-Matrix when the program is launched"
    with open(r'filter/FilterClonetypes.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    # initialising filter with database parameters
    # kmer20 = 115000000
    # kmer31 = 122000000
    BF = BF_v2.AbaumanniiBloomfilter(123000000)
    BF.set_arraysize(123000000)
    BF.set_hashes(7)
    BF.set_k(20)
    #paths = sorted(os.listdir(r"filter/species/"))
    paths = [r'filter/IC1.txt',
             r'filter/IC2.txt',
             r'filter/IC3.txt',
             r'filter/IC4.txt',
             r'filter/IC5.txt',
             r'filter/IC6.txt',
             r'filter/IC7.txt',
             r'filter/IC8.txt']
    BF.read_clonetypes(paths, clonetypes)
    return BF


def pre_processing():
    "Preprocesses the Bloomfilter-Matrix when the program is launched"
    with open(r'filter/FilterSpecies.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    # initialising filter with database parameters
    # kmer20 = 115000000
    # kmer31 = 122000000
    # kmer20+reversed 230000000
    # kmers21 173000000
    BF = BF_v2.AbaumanniiBloomfilter(115000000)
    BF.set_arraysize(115000000)
    BF.set_hashes(7)
    BF.set_k(21)
    #paths = sorted(os.listdir(r"filter/species_reversed/"))
    paths = sorted(os.listdir(r"filter/species/"))
    for i in range(len(paths)):
        paths[i] = r"filter/species/" + paths[i]
    BF.read_clonetypes(paths, clonetypes)
    return BF


# needs to be deleted i guess
def pre_processing_prefilter():
    "Preprocesses the Bloomfilter-Matrix when the program is launched"
    with open(r'filter/FilterHuman.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    # initialising filter with database parameters
    BF = BF_v2.AbaumanniiBloomfilter(23000000000)
    BF.set_arraysize(23000000000)
    BF.set_hashes(7)
    BF.set_k(20)
    paths = sorted(os.listdir(r"filter/human/"))
    for i in range(len(paths)):
        paths[i] = r"filter/human/" + paths[i]
    BF.read_clonetypes(paths, clonetypes)
    return BF


def pre_processing_prefilter2():
    "Preprocesses Acinetobacter Prefilter, collapse with other prefilter after testing"
    with open(r'filter/FilterAcinetobacter.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    # initialising filter with database parameters
    BF = BF_v2.AbaumanniiBloomfilter(3080000000)
    BF.set_arraysize(3080000000)
    BF.set_hashes(7)
    BF.set_k(21)
    paths = sorted(os.listdir(r"filter/acinetobacter/"))
    for i in range(len(paths)):
        paths[i] = r"filter/acinetobacter/" + paths[i]
    BF.read_clonetypes(paths, clonetypes)
    return BF


def pre_processing_prefilter_Culicidae():
    "Preprocesses the mosquito prefilter"
    with open(r'filter/FilterCulicidae.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    # initialising filter with database parameters
    # must be the same values used when trained!!
    BF = BF_v2.AbaumanniiBloomfilter(9000000)
    BF.set_arraysize(9000000)
    BF.set_hashes(7)
    BF.set_k(21)
    paths = sorted(os.listdir(r"filter/Culicidae_Prefilter/"))
    for i in range(len(paths)):
        paths[i] = r"filter/Culicidae_Prefilter/" + paths[i]
    BF.read_clonetypes(paths, clonetypes)
    return BF


def pre_processing_Culicidae_species():
    "Preprocesses the mosquito species filter"
    with open(r'filter/FilterCulicidaeSpecies.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    # initialising filter with database parameters
    # must be the same values used when trained!!
    BF = BF_v2.AbaumanniiBloomfilter(350000)
    BF.set_arraysize(350000)
    BF.set_hashes(7)
    BF.set_k(21)
    paths = sorted(os.listdir(r"filter/Culicidae_species/"))
    for i in range(len(paths)):
        paths[i] = r"filter/Culicidae_species/" + paths[i]
    BF.read_clonetypes(paths, clonetypes)
    return BF


def read_search_pre(reads, BF_pre, ext):
    reads_new = []
    counter = 0
    BF_pre.number_of_kmeres = 0
    BF_pre.hits_per_filter = [0]
    read_amount = 0
    reads_oxa_prefilter = []
    reads_oxa_filtered = []
    for single_read in reads:
        read_kmers = []
        hit_sum = sum(BF_pre.hits_per_filter)
        hits_per_filter_copy = BF_pre.hits_per_filter[:]
        # use a scaling sample size for contigs/scaffolds
        if ext == "fasta" or ext == "fna" or ext == "fa":
            sample_size = int(len(single_read) ** 0.5)
            threshold_read = sample_size * 0.7
            for i in range(0, len(single_read) - BF_pre.k, sample_size):
                if "N" not in single_read[i: i + BF_pre.k]:
                    BF_pre.lookup(single_read[i: i + BF_pre.k])
        # for reads use a static sample of 5
        # Taking sum of list as reference, if sum has not increased after testing those 3 kmeres,
        # then the read won't be tested further
        else:
            # TO-DO implement dynamic sample size
            k1 = single_read[0:BF_pre.k]  # first k-mer
            k2 = single_read[len(single_read) - BF_pre.k:]  # last k-mer
            mid = len(single_read) // 2
            k3 = single_read[mid:mid + BF_pre.k]  # k-mer in middle
            k4 = single_read[BF_pre.k:BF_pre.k * 2]
            k5 = single_read[mid + BF_pre.k:mid + BF_pre.k * 2]
            if "N" not in single_read:
                BF_pre.lookup(k1)
                BF_pre.lookup(k2)
                BF_pre.lookup(k3)
                BF_pre.lookup(k4)
                BF_pre.lookup(k5)
            threshold_read = 3
        # needs at least 2 of 3 hits to continue with read
        counter = 0
        if (sum(BF_pre.hits_per_filter) - hit_sum) > threshold_read:
            read_amount += 1
            for j in range(len(single_read) - BF_pre.k):
                if "N" not in single_read[j: j + BF_pre.k]:
                    read_kmers.append(single_read[j: j + BF_pre.k])
                    if ext == "fasta" or ext == "fna" or ext == "fa":
                        counter += 1
                        # extract up to 5000 kmeres per read/contig
                        if counter >= 5000:
                            break
            reads_oxa_prefilter.append(single_read)
            reads_new.append(read_kmers)
            BF_pre.hits_per_filter = hits_per_filter_copy
        else:
            # resetting hit counter
            BF_pre.hits_per_filter = hits_per_filter_copy
    reads_filtered = []
    threshold_dic = {}
    if ext == "fasta" or ext == "fna" or ext == "fa":
        cutoff = 0.7
    else:
        cutoff = 0.8
    counter = 0
    for i in range(len(reads_new)):
        threshold = 0
        for j in range(len(reads_new[i])):
            BF_pre.number_of_kmeres += 1
            hits_per_filter_copy = BF_pre.hits_per_filter[:]
            BF_pre.lookup(reads_new[i][j])
            if hits_per_filter_copy != BF_pre.hits_per_filter:
                threshold += 1
        if threshold >= cutoff * len(reads_new[i]):
            reads_filtered.append(reads_new[i])
            reads_oxa_filtered.append(reads_oxa_prefilter[i])
            counter += len(reads_new[i])
       # if ext == "fasta" or ext == "fna" or ext == "fa":
         #   if counter >= 50000:
           #     break
    return reads_filtered, reads_oxa_filtered


def read_search_spec(reads, quick, BF, ext):    
    "Searches sequence-data in Bloomfilter and gets kmer-hits"
    if quick < 4:
        BF.lookup_txt(reads, ext, quick)
        score = BF.get_score()
        hits = BF.get_hits_per_filter()
        names = BF.get_names()
        return score, names, hits, None
    # Metagenome mode
    elif quick == 4:
        reads_classified, predictions = BF.lookup_txt(reads, ext, quick)
        hits = None
        names = None
        return reads_classified, names, hits, predictions


def pre_processing_oxa():
    # getting filters
    oxa_families = sorted(os.listdir(r"filter/OXAs/families"))
    oxa_family_names = []
    paths = []
    for filter in oxa_families:
        oxa_family_names.append(filter[:-4])
    # get paths for oxa-family BF
    for i in range(len(oxa_families)):
        paths.append(r"filter/OXAs/families/" + oxa_families[i])

    # getting filters of individiual filters of oxa-families
    oxas_ind = sorted(os.listdir(r"filter/OXAs/individual"))
    # get paths for individiual BF
    paths_ind = []
    oxa_ind_names = []
    for i in range(len(oxas_ind)):
        paths_ind.append(r"filter/OXAs/individual/" + oxas_ind[i])
    # TODO: Rename variables
    paths_ind_ind = {}
    for i in range(len(paths_ind)):
        temp = sorted(os.listdir(paths_ind[i]))
        temp_list = []
        for j in range(len(temp)):
            temp_list.append(paths_ind[i] + "/" + temp[j])
        paths_ind_ind[oxas_ind[i]] = temp_list
    # list of BF-objects
    BF_dict = {}
    # initialising filter with database parameters
    BF = BF_v2.AbaumanniiBloomfilter(80000)
    BF.set_arraysize(80000)
    BF.set_clonetypes(len(paths))
    BF.set_hashes(7)
    BF.set_k(21)
    # User Options
    # reading single OXA filters
    BF.read_clonetypes(paths, oxa_family_names)
    BF_dict["OXA-families"] = BF
    # Add one BF-object for each oxa-family which contains individual oxa-BF
    for name, path_oxa_family in paths_ind_ind.items():
        names = []
        for filter in path_oxa_family:
            temp = filter.split("/")
            names.append(temp[-1][:-4])
        # initialising filter with database parameters
        BF = BF_v2.AbaumanniiBloomfilter(80000)
        BF.set_arraysize(80000)
        BF.set_clonetypes(len(path_oxa_family))
        BF.set_hashes(7)
        BF.set_k(21)
        # User Options
        # reading single OXA filters
        BF.read_clonetypes(path_oxa_family, names)
        BF_dict[name] = BF
    return BF_dict


def single_oxa(reads, ext, pipe=None):
    """Uses the Bloomfilter module to lookup the OXA-genes"""
    # getting filters
    paths = sorted(os.listdir(r"filter/OXAs/families/"))
    oxas = []
    for i in paths:
        oxas.append(i[:-4])

    for i in range(len(paths)):
        paths[i] = r"filter/OXAs/families/" + paths[i]


    # initialising filter with database parameters
    BF = BF_v2.AbaumanniiBloomfilter(80000)
    BF.set_arraysize(80000)
    BF.set_clonetypes(len(paths))
    BF.set_hashes(7)
    BF.set_k(21)
    # User Options

    # reading single OXA filters
    BF.read_clonetypes(paths, oxas)

    # starting Bloomfilter process, depends on filetype
    coordinates_forward, coordinates_reversed = BF.lookup_oxa(reads, ext)

    score = BF.get_oxa_score()
    BF.cleanup()
    del BF

    if pipe is not None:
        pipe.send([score, oxas])
        pipe.close()
    else:
        return score, oxas, coordinates_forward, coordinates_reversed


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
    start = time.time()
    p1 = Process(target=read_search, args=(IC_lookup, reads_ct, quick, child_ic))
    p1.start()
    p2 = Process(target=single_oxa, args=(reads, ext, child_oxa))
    p2.start()
    p1.join()
    p2.join()
    end = time.time()
    needed = round(end - start, 2)
    print("Time needed multiprocessing: ", needed)

    # getting results back from pipes
    results_ic = parent_ic.recv()  # has scores and names
    results_oxa = parent_oxa.recv()  # has scores and names

    return results_ic[0], results_ic[1], results_ic[2], results_oxa[0], results_oxa[1]
