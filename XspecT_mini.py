import search_filter
import os
from Bio import SeqIO, SeqRecord, Seq
from Bio.Seq import Seq
import Classifier
from OXA_Table import OXATable
import warnings
from copy import deepcopy
import time
import csv
import pickle
import matplotlib.pyplot as plt
import statistics
import random
from numpy import array
from numpy import sum
from collections import Counter
import psutil
import Bootstrap as bs
import sys
from pathlib import Path
#import map_kmers as map
from train_filter.interface_XspecT import load_translation_dict

from bitarray import bitarray

warnings.filterwarnings("ignore")


def xspecT_mini(file_path, XspecT, ClAssT, oxa, file_format, read_amount, csv_table, metagenome, genus, mode):
    """performs a BF-lookup for a set of genomes for testing purpose"""
    itemlist = ["albensis", "apis", "baretiae", "baumannii", "baylyi", "beijerinckii", "bereziniae",
                "bohemicus", "boissieri", "bouvetii", "brisouii", "calcoaceticus",
                "celticus", "chengduensis", "chinensis", "colistiniresistens", "courvalinii", "cumulans",
                "defluvii", "dispersus", "equi", "gandensis", "gerneri", "gs06", "gs16", "guerrae",
                "guillouiae", "gyllenbergii", "haemolyticus", "halotolerans", "harbinensis", "idrijaensis", "indicus",
                "johnsonii", "junii", "kanungonis", "kookii", "kyonggiensis", "lactucae", "lanii", "larvae",
                "lwoffii", "marinus", "modestus", "nectaris", "nosocomialis", "oleivorans", "parvus",
                "piscicola", "pittii", "pollinis", "populi", "portensis", "pseudolwoffii", "pullicarnis",
                "pragensis", "proteolyticus", "puyangensis",
                "qingfengensis", "radioresistens", "rathckeae", "rongchengensis", "rudis", "schindleri", "seifertii",
                "seohaensis", "shaoyimingii", "sichuanensis", "soli", "stercoris", "tandoii", "terrae",
                "terrestris", "tianfuensis", "tjernbergiae", "towneri", "ursingii", "variabilis", "venetianus",
                "vivanii", "wanghuae", "wuhouensis", "sp."]
    print("Preparing Bloomfilter...")
    start = time.time()
    if XspecT:
        #BF = search_filter.pre_processing()
        #Phillip
        # Getting the array sizes for pre processing of all bloomfilters.
        genera = search_filter.get_genera_array_sizes()

        # Pre processing of the bloomfilters for the species.
        BF = search_filter.pre_process_all(genera, k=21, meta_mode=False, genus=[genus])

         # aktuelle Speichernutzung auslesen
        process = psutil.Process()
        memory_info = process.memory_info()
        # Ausgabe des Speicherverbrauchs
        print(f"Aktueller Speicherverbrauch mit den Spezies BF: {memory_info.rss / 1024 / 1024:.2f} MB")

        # BF_1 = search_filter.pre_processing_prefilter()
        #BF_1_1 = search_filter.pre_processing_prefilter2()
        #Phillip
        # Pre processing of the bloomfilters for the metagenome mode.
        BF_1_1 = search_filter.pre_process_all(genera, k=21, meta_mode=True, genus=[genus])

         # aktuelle Speichernutzung auslesen
        process = psutil.Process()
        memory_info = process.memory_info()
        # Ausgabe des Speicherverbrauchs
        print(f"Aktueller Speicherverbrauch mit dem Master BF: {memory_info.rss / 1024 / 1024:.2f} MB")

    if ClAssT:
        BF_2 = search_filter.pre_processing_ClAssT()
    if oxa:
        BF_3 = search_filter.pre_processing_oxa()
    end = time.time()
    needed = round(end - start, 2)
    print("Time needed for preprocessing: ", needed)
    try:
        files = sorted(os.listdir(file_path))
    except FileNotFoundError:
        print("Error: Invalid filepath!")
        quit()
    if file_format == "fna" or file_format == "fasta" or file_format == "fa":
        for i in range(len(files) - 1, -1, -1):
            if 'fna' in files[i] or 'fasta' in files[i]:
                continue
            else:
                del files[i]
    elif file_format == "fastq" or file_format == "fq":
        for i in range(len(files) - 1, -1, -1):
            if 'fastq' in files[i] or 'fq' in files[i]:
                continue
            else:
                del files[i]
    if len(files) == 0:
        print("Error: No " + str(file_format) + " files in directory!")
        quit()
    paths = files[:]
    file_path2 = file_path[:]
    for i in range(len(file_path2)):
        if file_path2[i] == "\\":
            list_temp = list(file_path2)
            list_temp[i] = '/'
            file_path2 = ''.join(list_temp)
    start = time.time()
    for i in range(len(files)):
        paths[i] = file_path2 + "/" + paths[i]
    # Testing purpose delete later
    # extracts the GCF-Number
    # files_split = []
    # for i in range(len(files)):
    #    files_split.append(files[i].split("_"))
        # try:
    #        files_split[i] = files_split[i][0] + "_" + files_split[i][1]
    #    except:
    #        files_split[i] = files[i].split('.')[-2]
    # GCF_taxon = []
#    for i in range(len(files)):
#        with open(paths[i]) as file:
    #        head = file.readline()
    #        head = head.split()
    #        if head[2] == "sp.":
    #            GCF_taxon.append("none")
    #        else:
        #        GCF_taxon.append(head[2])
            # if head[2] != "baumannii":
                # del files[i]
                # del paths[i]
#    excelv3 = []
#    for i in range(len(files)):
    #    excelv3.append(files_split[i] + "," + GCF_taxon[i])
    # for i in range(0, len(excelv3)):
    #    excelv3[i] = [excelv3[i]]
#    with open(r'Results/XspecT_mini_csv/Test.csv', 'w', newline='') as file:
    #    writer = csv.writer(file)
    #    writer.writerows(excelv3)
    if XspecT:
        predictions, scores = xspecT(
            BF[genus], BF_1_1[genus], files, paths, file_format, read_amount, metagenome, genus, mode)
    if ClAssT:
        predictions_ClAssT, scores_ClAssT = clAssT(
            BF_2, files, paths, file_format, read_amount)
    if oxa:
        scores_oxa, scores_oxa_ind = blaOXA(BF_3, files, paths, file_format, read_amount)
    #if BioMonitoring:
        #predictions, scores = xspecT(
            #BF, BF_1_1, files, paths, file_format, read_amount, metagenome, BioMonitoring)
    print("Preparing results...")
    print("")
    end = time.time()
    needed = round(end - start, 2)
    print("Time needed: ", needed)
    print("")
    header_filename = "Filename"
    spaces = []
    space = "           "
    underscore = "________"
    name_max = len(max(itemlist, key=len))
    if XspecT:
        for i in range(len(predictions)):
            while len(predictions[i]) < name_max:
                predictions[i] += " "
    file_max = len(max(files, key=len))
    while len(header_filename) < file_max:
        header_filename += " "
        underscore += "_"
    for j in range(len(files)):
        for i in range(len(header_filename)-len(files[j])):
            space += " "
        spaces.append(space)
        space = "           "
    excel = []
    # formatting
    if ClAssT:
        for i in range(len(predictions_ClAssT)):
            if predictions_ClAssT[i] != "none" and predictions_ClAssT[i] != "None":
                predictions_ClAssT[i] += " "
    if XspecT and ClAssT:
        for i in range(len(scores_ClAssT)):
            if scores[i] == "1.0":
                scores[i] += " "

    if XspecT and ClAssT and oxa:
        excelv2 = []
        print(scores_oxa)
        print(scores_oxa_ind)
        for i in range(len(files)):
            if scores_oxa == ["None"]:
                excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i] + "           " +
                         predictions_ClAssT[i] + "            " + scores_ClAssT[i] + "           " + str(scores_oxa[i]) + "                             " + str(scores_oxa_ind[i][0]) +  "              " + str(scores_oxa_ind[i][1]))
            else:
                excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i] + "           " +
                            predictions_ClAssT[i] + "            " + scores_ClAssT[i] + "           " + str(scores_oxa[i]) + "           " + str(scores_oxa_ind[i][0]) +  "           " + str(scores_oxa_ind[i][1]))
            excelv2.append(files[i] + "," + predictions[i] + "," + scores[i] +
                           predictions_ClAssT[i] + "," + scores_ClAssT[i] + "," + str(scores_oxa[i]))
        print(header_filename + "           Species                  Score          Sub-Type        Score          blaOXA-Family                    blaOXA-Gene       Score")
        print(underscore +      "___________________________________________________________________________________________________________________________________________")
        for i in excel:
            print(i)
        for i in range(0, len(excelv2)):
            excelv2[i] = [excelv2[i]]
        if csv_table:
            with open(r'Results/XspecT_mini_csv/Results.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(excelv2)
        print("")
        print("")
    elif XspecT and not ClAssT and not oxa:
        excelv2 = []
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] +
                         predictions[i] + "       " + scores[i])
            excelv2.append(files[i] + "," + predictions[i] + "," + scores[i])
        print(header_filename + "           Species                  Score")
        print(underscore + "_________________________________________")
        for i in excel:
            print(i)
        for i in range(0, len(excelv2)):
            excelv2[i] = [excelv2[i]]
        if csv_table:
            with open(r'Results/XspecT_mini_csv/Results_XspecT.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(excelv2)
        print("")
        print("")
    elif ClAssT and not XspecT and not oxa:
        excelv2 = []
        for i in range(len(files)):
            excel.append(
                files[i] + spaces[i] + predictions_ClAssT[i] + "            " + scores_ClAssT[i])
            excelv2.append(
                files[i] + "," + predictions_ClAssT[i] + "," + scores_ClAssT[i])
        print(header_filename + "           Sub-Type        Score")
        print(underscore + "________________________________")
        for i in excel:
            print(i)
        print("")
        print("")
        for i in range(0, len(excelv2)):
            excelv2[i] = [excelv2[i]]
        if csv_table:
            with open(r'Results/XspecT_mini_csv/Results_ClAssT.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(excelv2)
    elif oxa and not ClAssT and not XspecT:
        excelv2 = []
        for i in range(len(files)):
            if scores_oxa == ["None"]:
                excel.append(files[i] + spaces[i] + str(scores_oxa[i]) + "                             " + str(scores_oxa_ind[i][0]) +  "              " + str(scores_oxa_ind[i][1]))
            else:
                excel.append(files[i] + spaces[i] + str(scores_oxa[i]) + "           " + str(scores_oxa_ind[i][0]) +  "           " + str(scores_oxa_ind[i][1]))

            excelv2.append(files[i] + "," + str(scores_oxa[i]))
        print(header_filename + "           blaOXA-Family                    blaOXA-Gene       Score")
        print(underscore + "_______________________________________________________________________")
        for i in excel:
            print(i)
        print("")
        print("")
        for i in range(0, len(excelv2)):
            excelv2[i] = [excelv2[i]]
        if csv_table:
            with open(r'Results/XspecT_mini_csv/Results_Oxa.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(excelv2)
    elif XspecT and ClAssT and not oxa:
        excelv2 = []
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i] +
                         "           " + predictions_ClAssT[i] + "            " + scores_ClAssT[i])
            excelv2.append(files[i] + "," + predictions[i] + "," + scores[i] +
                           "," + predictions_ClAssT[i] + "," + scores_ClAssT[i])
        print(header_filename +
              "           Species                  Score          Sub-Type        Score")
        print(underscore + "________________________________________________________________________")
        for i in excel:
            print(i)
        print("")
        print("")
        for i in range(0, len(excelv2)):
            excelv2[i] = [excelv2[i]]
        if csv_table:
            with open(r'Results/XspecT_mini_csv/Results.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(excelv2)
    elif XspecT and oxa and not ClAssT:
        excelv2 = []
        for i in range(len(files)):
            if scores_oxa == ["None"]:
                excel.append(files[i] + spaces[i] + predictions[i] +
                         "       " + scores[i] + "           " + str(scores_oxa[i]) + "                             " + str(scores_oxa_ind[i][0]) +  "              " + str(scores_oxa_ind[i][1]))
            else:
                excel.append(files[i] + spaces[i] + predictions[i] +
                         "       " + scores[i] + "           " + str(scores_oxa[i]) + "           " + str(scores_oxa_ind[i][0]) +  "           " + str(scores_oxa_ind[i][1]))
            excelv2.append(files[i] + "," + predictions[i] +
                           "," + scores[i] + str(scores_oxa[i]))
        print(header_filename +
              "           Species                  Score          blaOXA-Family                    blaOXA-Gene       Score")
        print(underscore + "_______________________________________________________________________________________________________________")
        for i in excel:
            print(i)
        for i in range(0, len(excelv2)):
            excelv2[i] = [excelv2[i]]
        if csv_table:
            with open(r'Results/XspecT_mini_csv/Results.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(excelv2)
        print("")
        print("")
    elif ClAssT and oxa and not XspecT:
        excelv2 = []
        for i in range(len(files)):
            if scores_oxa == ["None"]:
                excel.append(files[i] + spaces[i] + predictions_ClAssT[i] +
                         "            " + scores_ClAssT[i] + "           " + str(scores_oxa[i]) + "                             " + str(scores_oxa_ind[i][0]) +  "              " + str(scores_oxa_ind[i][1]))
            else:
                excel.append(files[i] + spaces[i] + predictions_ClAssT[i] +
                         "            " + scores_ClAssT[i] + "           " + str(scores_oxa[i]) + "           " + str(scores_oxa_ind[i][0]) +  "           " + str(scores_oxa_ind[i][1]))
            excelv2.append(files[i] + "," + predictions_ClAssT[i] +
                           "," + scores_ClAssT[i] + "," + str(scores_oxa[i]))
        print(header_filename +
              "           Sub-Type        Score          blaOXA-Family                    blaOXA-Gene       Score")
        print(underscore + "______________________________________________________________________________________________________")
        for i in excel:
            print(i)
        for i in range(0, len(excelv2)):
            excelv2[i] = [excelv2[i]]
        if csv_table:
            with open(r'Results/XspecT_mini_csv/Results.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(excelv2)
        print("")
        print("")


def xspecT(BF, BF_1_1, files, paths, file_format, read_amount, metagenome, genus, mode):
    """performs a BF-lookup for a set of genomes for testing purpose"""
    matrix: bitarray
    matrix = BF_1_1.matrix
    ones = matrix.count(1)
    zeroes = matrix.count(0)
    print(f"1: {ones}, 0: {zeroes}")
    print(f"len: {len(matrix)}, sum: {ones+zeroes}")

    print("Starting taxonomic assignment on species-level...")
    predictions = []
    scores = []
    counterx = 0
    contig_header = []
    contig_seq = []
    # Phillip
    names_path = Path(__file__).parent.absolute() / 'filter' / 'species_names' / ('Filter' + genus + '.txt')
    with open(names_path, 'rb') as fp:
        names = pickle.load(fp)
    names = sorted(names)
    translation_dict = load_translation_dict(genus)
    #
    for i in range(len(files)):
        if i == int(len(files)/6) or i == int(len(files)/3) or i == int(len(files)/2) or i == int(len(files)/1.5) or i == int(len(files)/1.2):
            print("...")
        BF.number_of_kmeres = 0
        BF.hits_per_filter = [0] * BF.clonetypes
        BF_1_1.number_of_kmeres = 0
        BF_1_1.hits_per_filter = [0]
        if file_format == "fasta" or file_format == "fna" or file_format =="fa":
            if metagenome:               
                contigs = []
                contigs_classified = {}
                for sequence in SeqIO.parse(paths[i], "fasta"):
                    contigs = []
                    contigs_kmers = []
                    BF_1_1.kmer_hits_single = []
                    BF_1_1.number_of_kmeres = 0
                    BF_1_1.hits_per_filter = [0] * BF.clonetypes
                    # Taking sum of list as reference, if sum has not increased after testing those 3 kmeres,
                    # then the contigs won't be tested further
                    hit_sum = sum(BF_1_1.hits_per_filter)
                    hits_per_filter_copy = BF_1_1.hits_per_filter[:]
                    sample_size = int(len(str(sequence.seq)) ** 0.5)
                    threshold_contig = sample_size * 0.7
                    for i in range(0, len(str(sequence.seq)) - BF_1_1.k, sample_size):
                        if "N" not in str(sequence.seq[i: i + BF_1_1.k]):
                            BF_1_1.lookup(str(sequence.seq[i: i + BF_1_1.k]).upper())

		            # needs at least 70% hits to continue with the contig
                    counter = 0
                    print(f"Sum: {sum(BF_1_1.hits_per_filter)}; threshold: {threshold_contig}")
                    if (sum(BF_1_1.hits_per_filter) - hit_sum) > threshold_contig:
                        for j in range(len(str(sequence.seq)) - BF_1_1.k):
                            if "N" not in str(sequence.seq[j: j + BF_1_1.k]):
                                contigs_kmers.append(str(sequence.seq[j: j + BF_1_1.k]).upper())
                                counter += 1
                                # how many kmers? to use
                                if counter >= 5000000:
                                    break
		                # contigs_kmers.append(str(reverse_sequence[j: j + BF_1_1.k]))
                        contigs.append(contigs_kmers)
                        BF_1_1.hits_per_filter = hits_per_filter_copy
                    else:
		                # resetting hit counter
                        BF_1_1.hits_per_filter = hits_per_filter_copy
                        continue

                    contigs_filtered = []
                    counter = 0
                    # Since we classify individual contigs now, the var contigs only contains one item which makes those loops unneccesary
                    for i in range(len(contigs)):
                        threshold = 0
                        for j in range(len(contigs[i])):
                            BF_1_1.number_of_kmeres += 1
                            hits_per_filter_copy = BF_1_1.hits_per_filter[:]
                            BF_1_1.lookup(contigs[i][j])
                            if hits_per_filter_copy != BF_1_1.hits_per_filter:
                                threshold += 1
                        # parameter value needs to be determined
                        if threshold >= (0.7*len(contigs[i])):
                            contigs_filtered += contigs[i]
                            counter += len(contigs[i])
                        if counter >= 5000:
                            break

                    # since we do indv. contig classifications we need to reset the BF vars
                    BF.kmer_hits_single = []
                    BF.number_of_kmeres = 0
                    BF.hits_per_filter = [0] * BF.clonetypes
                    for kmer in contigs_filtered:
                        BF.number_of_kmeres += 1
                        kmer_reversed = str(Seq(kmer).reverse_complement())
                        if kmer > kmer_reversed:
                            BF.lookup(kmer)
                        else:
                            BF.lookup(kmer_reversed)
                    score = BF.get_score()
                    score_edit = [str(x) for x in score]
                    score_edit = ",".join(score_edit)

                    # making prediction
                    index_result = max(range(len(score)), key=score.__getitem__)
                    prediction = names[index_result]

                    # skip ambiguous contigs
                    if max(score) == sorted(score)[-2]:
                        continue

                    # bootstrapping
                    bootstrap_n = 100
                    samples = bs.bootstrap(BF.kmer_hits_single, BF.number_of_kmeres, bootstrap_n)
                    sample_scores = bs.bootstrap_scores(samples, BF.number_of_kmeres, BF.clonetypes)
                    bootstrap_score = 0
                    bootstrap_predictions = []
                    for i in range(len(sample_scores)):
                        # skip ambiguous contigs (species with same score)
                        if max(sample_scores[i]) != sorted(sample_scores[i])[-2]:
                            bootstrap_predictions.append(names[max(range(len(sample_scores[i])), key=sample_scores[i].__getitem__)])
                            if max(range(len(sample_scores[i])), key=sample_scores[i].__getitem__) == index_result:
                                bootstrap_score += 1
                        else:
                            continue
                    bootstrap_score = bootstrap_score/bootstrap_n

                    # ---------------------------------------------------------------------------------------------
                    # Collect results
                    # change this var to the species you want your contigs saved from
                    save_contigs = "none"

                    if (genus[0] + ". " + prediction) not in contigs_classified:
                        contigs_classified[genus[0] + ". " + prediction] = [[max(score)], 1, [len(str(sequence.seq))], sorted(score)[-2]/max(score), [bootstrap_score], contigs_filtered, None]
                        if prediction == save_contigs:
                            contig_header += [sequence.description]
                            contig_seq += [str(sequence.seq)]
                    else:
                        contigs_classified[genus[0] + ". " + prediction][0] += [max(score)]
                        contigs_classified[genus[0] + ". " + prediction][1] += 1
                        contigs_classified[genus[0] + ". " + prediction][2] += [len(str(sequence.seq))]
                        contigs_classified[genus[0] + ". " + prediction][3] += sorted(score)[-2]/max(score)
                        contigs_classified[genus[0] + ". " + prediction][4] += [bootstrap_score]
                        contigs_classified[genus[0] + ". " + prediction][5] += contigs_filtered
                        if prediction == save_contigs:
                            contig_header += [sequence.description]
                            contig_seq += [str(sequence.seq)]
                        #scores.append(str(max(score)))
            else:
                # Important! Resetting the kmer_hits_single otherwise MEMORY LEAK
                BF.kmer_hits_single = []
                for sequence in SeqIO.parse(paths[i], "fasta"):
                    for j in range(0, len(sequence.seq) - BF.k, mode):
                        BF.number_of_kmeres += 1
                        kmer = str(sequence.seq[j : j + BF.k])
                        kmer_reversed = str(Seq(kmer).reverse_complement())
                        if kmer > kmer_reversed:
                            BF.lookup(kmer)
                        else:
                            BF.lookup(kmer_reversed)

            score = BF.get_score()
            #print("Scores: ", score)
            if metagenome:
                # map kmers to genome for HGT detection
                # change later to new functions this is OLD
                if False:
                    for prediction in contigs_classified:
                        kmers = contigs_classified[prediction][5]
                        # Strip "A."
                        prediction = prediction[2:]
                        # kmer mapping to genome, start by loading the kmer_dict in
                        path_pos = "filter\kmer_positions\Acinetobacter\\" + prediction + "_positions.txt"
                        # delete later
                        path_posv2 = "filter\kmer_positions\Acinetobacter\\" + prediction + "_complete_positions.txt"
                        # cluster kmers to contigs
                        # delete try later
                        try:
                            with open(path_pos, 'rb') as fp:
                                kmer_dict = pickle.load(fp)
                        except:
                            with open(path_posv2, 'rb') as fp:
                                kmer_dict = pickle.load(fp)
                        contig_amounts_distances = bs.cluster_kmers(kmers, kmer_dict)
                        contigs_classified[genus[0] + ". " + prediction][6] = contig_amounts_distances
                        #del kmer_dict
                for key, value in contigs_classified.items():
                    number_of_contigs = value[1]
                     # save results
                    results_clustering = [[key + "," + str(statistics.median(value[0])) + "," + str(number_of_contigs), str(statistics.median(value[2])) + "," + str(round(value[3]/number_of_contigs, 2)) + "," + str(statistics.median(value[4])) + "," + str(value[6])]]
                    #with open(r'Results/XspecT_mini_csv/Results_Clustering.csv', 'a', newline='') as file:
                        #writer = csv.writer(file)
                        #writer.writerows(results_clustering)
                    value[0] = "Score Median: " + str(statistics.median(value[0]))
                    value[1] = "Number of Contigs: " + str(number_of_contigs)
                    value[2] = "Contig-Length Median: " + str(statistics.median(value[2]))
                    value[3] = "Repetiviness: " + str(round(value[3]/number_of_contigs, 2))
                    value[4] = "Bootstrap Median: " + str(statistics.median(value[4]))
                    #value[6] = "Clusters: " + str(value[6])
                    contigs_classified[key] = value
                    print("Species: ", key)
                    print(value[0])
                    print(value[1])
                    print(value[2])
                    print(value[3])
                    print(value[4])
                    print(value[6])
                    print()

                save_contigs = "none"
                if save_contigs != "none":
                    with open(r"Results/Contigs_saved.fasta", 'w') as file:
                        for j in range(len(contig_header)):
                            file.write(contig_header[j] + "\n")
                            file.write(contig_seq[j] + "\n")
                            file.write("\n")
        elif file_format == "fastq" or file_format == "fq":
            if metagenome:
                # ---------------------------------------------------------------------------------------------
                # initialize variables
                BF_1_1.kmer_hits_single = []
                BF_1_1.number_of_kmeres = 0
                BF_1_1.hits_per_filter = [0] * BF.clonetypes
                counter = 0
                reads = []
                reads_classified = {}
                reads_passed = 0
                ambiguous_reads = 0

                # ---------------------------------------------------------------------------------------------
                # First prefiltering step: Check if read contains at least 3 kmeres
                for sequence in SeqIO.parse(paths[i], "fastq"):
                    dna_composition = {}
                    dna_composition = calculate_dna_composition(sequence.seq)
                    BF_1_1.kmer_hits_single = []
                    BF_1_1.number_of_kmeres = 0
                    BF_1_1.hits_per_filter = [0] * BF.clonetypes
		            # reverse_sequence = sequence.seq.reverse_complement()
                    read_kmers = []
                    reads = []
                    if counter < read_amount:
                        counter += 1
                    else:
                        break
                    k1 = str(sequence.seq[0:BF_1_1.k])  # first k-mer
                    k2 = str(sequence.seq[len(str(sequence.seq)) - BF_1_1.k:])  # last k-mer
                    mid = len(str(sequence.seq)) // 2
                    k3 = str(sequence.seq[mid:mid + BF_1_1.k])  # k-mer in middle
                    k4 = str(sequence.seq[BF_1_1.k:BF_1_1.k * 2])
                    k5 = str(sequence.seq[mid + BF_1_1.k:mid + BF_1_1.k * 2])
                    # Taking sum of list as reference, if sum has not increased after testing those 3 kmeres,
                    # then the read won't be tested further
                    hit_sum = sum(BF_1_1.hits_per_filter)
                    hits_per_filter_copy = BF_1_1.hits_per_filter[:]
                    #sample_size = int(len(str(sequence.seq)) ** 0.5)
                    #threshold_read = sample_size * 0.7
                    #for i in range(0, len(str(sequence.seq)) - BF_1_1.k, sample_size):
                    #    if "N" not in str(sequence.seq[i: i + BF_1_1.k]):
                    #        BF_1_1.lookup(str(sequence.seq[i: i + BF_1_1.k]))
                    if "N" not in str(sequence.seq):
                        BF_1_1.lookup(k1)
                        BF_1_1.lookup(k2)
                        BF_1_1.lookup(k3)
                        BF_1_1.lookup(k4)
                        BF_1_1.lookup(k5)
                    else:
                        continue
		            # needs at least 2 of 3 hits to continue with read
                    if (sum(BF_1_1.hits_per_filter) - hit_sum) > 3:
                        for j in range(len(str(sequence.seq)) - BF_1_1.k):
                            read_kmers.append(str(sequence.seq[j: j + BF_1_1.k]))
		                # read_kmers.append(str(reverse_sequence[j: j + BF_1_1.k]))
                        reads.append(read_kmers)
                        BF_1_1.hits_per_filter = hits_per_filter_copy
                    else:
		                # resetting hit counter
                        BF_1_1.hits_per_filter = hits_per_filter_copy
                        continue

                    # ---------------------------------------------------------------------------------------------
                    # Second prefiltering step: Check if read contains at least 80% of kmers from one species
                    #reads_filtered = set()
                    reads_filtered = []
                    for i in range(len(reads)):
                        threshold = 0
                        for j in range(len(reads[i])):
                            BF_1_1.number_of_kmeres += 1
                            hits_per_filter_copy = BF_1_1.hits_per_filter[:]
                            if "N" not in reads[i][j]:
                                BF_1_1.lookup(reads[i][j])
                            if hits_per_filter_copy != BF_1_1.hits_per_filter:
                                threshold += 1
                        if threshold >= 0.7 * len(reads[i]):
                            reads_filtered += reads[i]
                    if len(reads_filtered) == 0:
                        continue
                    
                    # ---------------------------------------------------------------------------------------------
                    # Start of the actual classification 
                    BF.number_of_kmeres = 0
                    BF.hits_per_filter = [0] * BF.clonetypes
                    BF.kmer_hits_single = []
                    for kmer in reads_filtered:
                        if "N" not in kmer:
                            BF.number_of_kmeres += 1
                            kmer_reversed = str(Seq(kmer).reverse_complement())
                            if kmer > kmer_reversed:
                                BF.lookup(kmer)
                            else:
                                BF.lookup(kmer_reversed)
                        else:
                            continue
                    score = BF.get_score()
                    score_edit = [str(x) for x in score]
                    score_edit = ",".join(score_edit)

                    # making prediction
                    index_result = max(range(len(score)), key=score.__getitem__)
                    prediction = names[index_result]
                    if max(score) == sorted(score)[-2]:
                        ambiguous_reads += 1
                        #print("Ambiguous read")
                        #continue

                    # ---------------------------------------------------------------------------------------------
                    # bootstrapping
                    bootstrap_n = 100
                    samples = bs.bootstrap(BF.kmer_hits_single, BF.number_of_kmeres, bootstrap_n)
                    sample_scores = bs.bootstrap_scores(samples, BF.number_of_kmeres, BF.clonetypes)
                    bootstrap_score = 0
                    bootstrap_predictions = []
                    for i in range(len(sample_scores)):
                        if max(sample_scores[i]) != sorted(sample_scores[i])[-2]:
                            bootstrap_predictions.append(names[max(range(len(sample_scores[i])), key=sample_scores[i].__getitem__)])
                            if max(range(len(sample_scores[i])), key=sample_scores[i].__getitem__) == index_result:
                                bootstrap_score += 1
                        else:
                            continue
                    bootstrap_score = bootstrap_score/bootstrap_n

                    # ---------------------------------------------------------------------------------------------
                    # HGT identification pipeline start


                    # identify split reads from HGT
                    #split_regions = map.identify_split_reads(score, BF.kmer_hits_single)

                    # split_read contains touples --> ([start, end], index of species)
                    
                    # get involed species
                    #split_reads_species = []
                    #if split_regions:
                    #    split_reads_species.append(names[split_regions[0][1]])
                    #    split_reads_species.append(names[split_regions[1][1]])


                    # ---------------------------------------------------------------------------------------------
                    # Collect results from classification
                    if (genus[0] + ". " + prediction) not in reads_classified:
                        reads_classified[genus[0] + ". " + prediction] = [max(score), 1, sorted(score)[-2]/max(score),  BF.number_of_kmeres, [bootstrap_score], reads_filtered, None, [dna_composition]]
                    else:
                        reads_classified[genus[0] + ". " + prediction][1] += 1
                        reads_classified[genus[0] + ". " + prediction][0] += max(score)
                        reads_classified[genus[0] + ". " + prediction][2] += sorted(score)[-2]/max(score)
                        reads_classified[genus[0] + ". " + prediction][3] += BF.number_of_kmeres
                        reads_classified[genus[0] + ". " + prediction][4] += [bootstrap_score]
                        reads_classified[genus[0] + ". " + prediction][5] += reads_filtered
                        reads_classified[genus[0] + ". " + prediction][7] += [dna_composition]

            else:
                # classification for sequence pure reads, check every 10th kmer (or everyone for "complete" mode)
                counter = 0
                # Important! Resetting the kmer_hits_single otherwise MEMORY LEAK
                BF.kmer_hits_single = []
                for sequence in SeqIO.parse(paths[i], "fastq"):
                    if counter < read_amount:
                        counter += 1
                        for j in range(0, len(sequence.seq) - BF.k+1, mode):
                            BF.number_of_kmeres += 1
                            kmer = str(sequence.seq[j : j + BF.k])
                            kmer_reversed = str(Seq(kmer).reverse_complement())
                            if kmer > kmer_reversed:
                                BF.lookup(kmer)
                            else:
                                BF.lookup(kmer_reversed)
                    else:
                        break
            score = BF.get_score()

            # ---------------------------------------------------------------------------------------------
            # Collect results from classification
            if metagenome:
                for key, value in reads_classified.items():
                    if key == "unknown":
                        continue
                    value.insert(2, value[0]/value[1])
                    value.pop(0)
                    reads_classified[key] = value
                    print(key, value[0], round(value[1], 2), round(value[2]/value[0], 2), round(value[3]/value[0], 2), statistics.median(value[4]), value[6], value[7])
            score_edit = [str(x) for x in score]
            score_edit = ",".join(score_edit)
        # making prediction
        if not metagenome:
            #prediction = Classifier.classify(r'Training_data/Training_data_spec.csv', score, True)
            # Phillip
            file_name = genus + "_Training_data_spec.csv"
            path = Path(__file__).parent.absolute() / "Training_data" / file_name
            prediction = Classifier.classify(path, score, True)
            # SVM TURNED OFF TEsting!!
            index_result = max(range(len(score)), key=score.__getitem__)
            #prediction = names[index_result]
            names_copy = names[:]
            # sort score by descending order and names_copy accordingly
            score, names_copy = zip(*sorted(zip(score, names_copy), reverse=True))
            #print(score[0:3])
            #print(names_copy[0:3])
        else:
            index_result = max(range(len(score)), key=score.__getitem__)
            prediction = names[index_result]
        # Phillip
        unknown_predictions = ["sp.", "unknown"]
        if prediction not in unknown_predictions:
            prediction_name = translation_dict[prediction]
            short_prediction_name = " ".join(prediction_name.split(" ")[1:])
            predictions.append(f"{genus[0]}. {short_prediction_name}")
        else:
            predictions.append(f"{genus[0]}. {prediction}")
        #
        scores.append(str(max(score)))
    print("Taxonomic assignment done...")
    return predictions, scores


def clAssT(BF_2, files, paths, file_format, read_amount):
    print("Starting strain-typing on sub-type-level...")
    predictions_ClAssT = []
    scores_ClAssT = []
    for i in range(len(files)):
        if i == int(len(files)/6) or i == int(len(files)/3) or i == int(len(files)/2) or i == int(len(files)/1.5) or i == int(len(files)/1.2):
            print("...")
        BF_2.number_of_kmeres = 0
        BF_2.hits_per_filter = [0] * BF_2.clonetypes
        if file_format == "fasta" or file_format == "fna":
            for sequence in SeqIO.parse(paths[i], "fasta"):
                # Originally 10
                for j in range(0, len(sequence.seq) - BF_2.k, 500):
                    BF_2.number_of_kmeres += 1
                    BF_2.lookup(str(sequence.seq[j: j + BF_2.k]))
        elif file_format == "fastq" or file_format == "fq":
            counter = 0
            for sequence in SeqIO.parse(paths[i], "fastq"):
                if counter < read_amount:
                    counter += 1
                    for j in range(0, len(sequence.seq) - BF_2.k+1, 10):
                        BF_2.number_of_kmeres += 1
                        BF_2.lookup(str(sequence.seq[j: j + BF_2.k]))
                else:
                    break
        score_ClAssT = BF_2.get_score()
        score_edit_ClAssT = [str(x) for x in score_ClAssT]
        score_edit_ClAssT = ",".join(score_edit_ClAssT)
        prediction_ClAssT = Classifier.classify(r'Training_data/Training_data_IC.csv', score_ClAssT, [True,True,True,True,True,True,True,True,False])
        predictions_ClAssT.append(prediction_ClAssT)
        scores_ClAssT.append(str(max(score_ClAssT)))

    print("Strain-typing on sub-type-level done...")
    return predictions_ClAssT, scores_ClAssT


def blaOXA(BF_3, files, paths, file_format, read_amount):
    start = time.time()
    print("Start screening for blaOXA-genes...")
    paths_oxa = sorted(os.listdir(r"filter/OXAs/families"))
    BF_families = BF_3["OXA-families"]
    oxas = []
    scores_oxa = []
    scores_oxa_ind = []
    for i in paths_oxa:
        oxas.append(i[:-4])
    #print("OXA-families: ", oxas) # correct
    for i in range(len(files)):
        oxa_dic = {}
        if i == int(len(files)/6) or i == int(len(files)/3) or i == int(len(files)/2) or i == int(len(files)/1.5) or i == int(len(files)/1.2):
            print("...")
        # Checking file type
        # if the file is fasta -> concat lines
        reads = []
        BF_families.number_of_kmeres = 0
        BF_families.hits_per_filter = [0] * BF_families.clonetypes
        BF_families.table = OXATable()
        BF_families.table.read_dic(r'filter/OXAs_dict/oxa_dict.txt')
        if file_format == "fasta" or file_format == "fna":
            for sequence in SeqIO.parse(paths[i], "fasta"):
                reads.append(str(sequence.seq))
            BF_families.lookup_oxa(reads, ".fna")
        elif file_format == "fastq" or file_format == "fq":
            counter = 0
            for sequence in SeqIO.parse(paths[i], "fastq"):
                if counter < read_amount:
                    counter += 1
                    reads.append(str(sequence.seq))
                else:
                    break
            BF_families.lookup_oxa(reads, ".fq")
            #print("Reads used: ", counter)
        score_oxa = BF_families.get_oxa_score()
        #print("Score: ", score_oxa)
        for i in range(len(oxas)):
            oxa_dic[oxas[i]] = score_oxa[i]
        for i in range(len(oxa_dic)):
            if oxa_dic[oxas[i]] < 0.3:
                del oxa_dic[oxas[i]]
        if len(oxa_dic) == 0:
            oxa_dic = "None"
        if oxa_dic != "None":
            oxa_dic = dict(sorted(oxa_dic.items(), key=lambda item: item[1]))
        scores_oxa.append(oxa_dic)
        # prepare data for next taxonomic level
        oxa_names = []
        #print(oxa_dic)
        for oxa_family in oxa_dic:
            oxa_names.append(oxa_family[:-7])
        for oxa_family in oxa_names:
            if oxa_dic == "None":
                scores_oxa_ind.append(["None", 0])
                break
            #print("blaOXA: ", oxa_dic)
            oxa_dic_ind = {}
            ## TODO:
            #print("blaOXA: ", oxa_family)
            BF_ind = BF_3[oxa_family]
            BF_ind.number_of_kmeres = 0
            BF_ind.hits_per_filter = [0] * BF_ind.clonetypes
            BF_ind.table = OXATable()
            BF_ind.table.read_dic(r'filter/OXAs_dict/oxa_dict.txt')
            paths_oxa = sorted(os.listdir(r"filter/OXAs/individual/" + oxa_family))
            oxas_ind = []
            for i in paths_oxa:
                oxas_ind.append(i[:-4])
            if file_format == "fasta" or file_format == "fna":
                BF_ind.lookup_oxa(reads, ".fna")
            elif file_format == "fastq" or file_format == "fq":
                BF_ind.lookup_oxa(reads, ".fq")
            score_oxa = BF_ind.get_oxa_score()
            # build dict with oxa-gen and its score
            for i in range(len(oxas_ind)):
                oxa_dic_ind[oxas_ind[i]] = score_oxa[i]
            # filter dict by score
            if len(oxa_dic_ind) == 0 or max(oxa_dic_ind.values()) < 0.3:
                scores_oxa_ind.append("None")
            else:
                scores_oxa_ind.append([max(oxa_dic_ind, key=oxa_dic_ind.get), oxa_dic_ind[max(oxa_dic_ind, key=oxa_dic_ind.get)]])
    end = time.time()
    needed = round(end - start, 2)
    print("Time needed: ", needed)
    print("Screening for blaOXA-genes done...")
    return scores_oxa, scores_oxa_ind


# Funktion zur Bestimmung der DNA-Zusammensetzung
def calculate_dna_composition(sequence):
    composition = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    total = 0
    for base in sequence:
        if base in composition:
            composition[base] += 1
            total += 1
    for base in composition:
        composition[base] = round(composition[base] / total, 2)
    
    return composition


def main():

    arg_list = sys.argv
    #Phillip
    genus = arg_list[1]
    genera = search_filter.get_genera_array_sizes()
    genera = list(genera.keys())

    # check for args that the tool doesnt know
    # list of all known args
    known_args = ["XspecT", "xspect", "ClAssT", "classt", "Oxa", "oxa", "Metagenome", "metagenome", "fasta", "fna", "fa", "fastq", "fq", "save", "Save", "complete", "Complete"]
    # TODO

    if genus not in genera:
        print(f"{genus} is unknown.")
        quit()
    if "XspecT" in arg_list or "xspect" in arg_list:
        XspecT = True
    else:
        XspecT = False
    if "ClAssT" in arg_list or "classt" in arg_list and genus == "Acinetobacter":
        ClAssT = True
    elif "ClAssT" in arg_list or "classt" in arg_list and genus != "Acinetobacter":
        print(f"ClAssT unavailable for {genus}")
    else:
        ClAssT = False
    if "Oxa" in arg_list or "oxa" in arg_list and genus == "Acinetobacter":
        oxa = True
    elif "Oxa" in arg_list or "oxa" in arg_list and genus != "Acinetobacter":
        print(f"Oxa unavailable for {genus}")
    else:
        oxa = False
    if "Metagenome" in arg_list or "metagenome" in arg_list:
        metagenome = True
    else:
        metagenome = False
    if ("fasta" in arg_list) or ("fna" in arg_list) or ("fa" in arg_list):
        file_format = "fasta"
        read_amount = 342480
    elif ("fastq" in arg_list) or ("fq" in arg_list):
        file_format = "fastq"
        index = arg_list.index("fastq")
        if arg_list[index+1].isdigit():
            read_amount = int(arg_list[index+1])
        else:
            print("Error: Wrong Input, use a number after fastq!")
            quit()
    else:
        print("Error: Wrong Input, use fasta/fna/fa or fastq/fq!")
        quit()
    if "save" in arg_list or "Save" in arg_list:
        csv_table = True
    else:
        csv_table = False
    if "complete" in arg_list or "Complete" in arg_list:
        mode = 1
    else:
        mode = 500

    file_path = arg_list[-1]
    
    xspecT_mini(file_path, XspecT, ClAssT, oxa, file_format, read_amount, csv_table, metagenome, genus, mode)


if __name__ == '__main__':
    main()
