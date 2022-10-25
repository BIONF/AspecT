import search_filter
import os
from Bio import SeqIO, SeqRecord, Seq
import Classifier
from OXA_Table import OXATable
import warnings
from copy import deepcopy
import time
import csv
import pickle
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

def xspecT_mini(file_path, XspecT, ClAssT, oxa, file_format, read_amount, csv_table, metagenome):
    """performs a BF-lookup for a set of genomes for testing purpose"""
    itemlist = ["albensis", "apis", "baretiae", "baumannii", "baylyi", "beijerinckii", "bereziniae",
                "bohemicus", "boissieri", "bouvetii", "brisouii", "calcoaceticus",
                "celticus", "chengduensis", "chinensis", "colistiniresistens","courvalinii", "cumulans",
                "defluvii", "dispersus", "equi", "gandensis", "gerneri","gs06","gs16", "guerrae",
                "guillouiae", "gyllenbergii", "haemolyticus", "halotolerans", "harbinensis", "idrijaensis", "indicus",
                "johnsonii", "junii", "kanungonis", "kookii", "kyonggiensis", "lactucae", "lanii", "larvae",
                "lwoffii", "marinus", "modestus", "nectaris", "nosocomialis", "oleivorans", "parvus",
                "piscicola", "pittii", "pollinis", "populi", "portensis", "pseudolwoffii", "pullicarnis",
                "pragensis", "proteolyticus","puyangensis",
                "qingfengensis", "radioresistens", "rathckeae", "rongchengensis", "rudis", "schindleri", "seifertii",
                "seohaensis", "shaoyimingii", "sichuanensis", "soli", "stercoris", "tandoii", "terrae",
                "terrestris", "tianfuensis", "tjernbergiae", "towneri", "ursingii","variabilis", "venetianus",
                "vivianii", "wanghuae", "wuhouensis", "sp."]
    print("Preparing Bloomfilter...")
    start = time.time()
    if XspecT:
        BF = search_filter.pre_processing()
        #BF_1 = search_filter.pre_processing_prefilter()
        BF_1_1 = search_filter.pre_processing_prefilter2()
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
    if file_format == "fna" or file_format =="fasta":
        for i in range(len(files) -1, -1, -1):
            if 'fna' in files[i] or 'fasta' in files[i]:
                continue
            else:
                del files[i]
    elif file_format == "fastq" or file_format =="fq":
        for i in range(len(files) -1, -1, -1):
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
    #extracts the GCF-Number
    #files_split = []
    #for i in range(len(files)):
    #    files_split.append(files[i].split("_"))
        #try:
    #        files_split[i] = files_split[i][0] + "_" + files_split[i][1]
    #    except:
    #        files_split[i] = files[i].split('.')[-2]
    #GCF_taxon = []
#    for i in range(len(files)):
#        with open(paths[i]) as file:
    #        head = file.readline()
    #        head = head.split()
    #        if head[2] == "sp.":
    #            GCF_taxon.append("none")
    #        else:
        #        GCF_taxon.append(head[2])
            #if head[2] != "baumannii":
                #del files[i]
                #del paths[i]
#    excelv3 = []
#    for i in range(len(files)):
    #    excelv3.append(files_split[i] + "," + GCF_taxon[i])
    #for i in range(0, len(excelv3)):
    #    excelv3[i] = [excelv3[i]]
#    with open(r'Results/XspecT_mini_csv/Test.csv', 'w', newline='') as file:
    #    writer = csv.writer(file)
    #    writer.writerows(excelv3)
    if XspecT:
        predictions, scores = xspecT(BF, BF_1_1, files, paths, file_format, read_amount, metagenome)
    if ClAssT:
        predictions_ClAssT, scores_ClAssT = clAssT(BF_2, files, paths, file_format, read_amount)
    if oxa:
        scores_oxa = blaOXA(BF_3, files, paths, file_format, read_amount)
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
    #formatting
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
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i] + "           " + predictions_ClAssT[i] + "            " + scores_ClAssT[i] + "           " + str(scores_oxa[i]))
            excelv2.append(files[i] + "," + predictions[i] + "," + scores[i] + predictions_ClAssT[i] + "," + scores_ClAssT[i] + "," + str(scores_oxa[i]))
        print(header_filename + "           Species                  Score          Sub-Type        Score          blaOXA-Genes")
        print(underscore + "__________________________________________________________________________________________________")
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
            excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i])
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
            excel.append(files[i] + spaces[i] + predictions_ClAssT[i] + "            " + scores_ClAssT[i])
            excelv2.append(files[i] + "," + predictions_ClAssT[i] + "," + scores_ClAssT[i])
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
            excel.append(files[i] + spaces[i] + str(scores_oxa[i]))
            excelv2.append(files[i] + "," + str(scores_oxa[i]))
        print(header_filename + "           blaOXA-Genes")
        print(underscore + "___________________________")
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
            excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i] + "           " + predictions_ClAssT[i] + "            " + scores_ClAssT[i])
            excelv2.append(files[i] + "," + predictions[i] + "," + scores[i] + "," + predictions_ClAssT[i] + "," + scores_ClAssT[i])
        print(header_filename + "           Species                  Score          Sub-Type        Score")
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
            excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i]  + "           " + str(scores_oxa[i]))
            excelv2.append(files[i] + "," + predictions[i] + "," + scores[i] + str(scores_oxa[i]))
        print(header_filename + "           Species                  Score          blaOXA-Genes")
        print(underscore + "___________________________________________________________________")
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
            excel.append(files[i] + spaces[i] + predictions_ClAssT[i] + "            " + scores_ClAssT[i] + "           " + str(scores_oxa[i]))
            excelv2.append(files[i] + "," + predictions_ClAssT[i] + "," + scores_ClAssT[i] + "," + str(scores_oxa[i]))
        print(header_filename + "           Sub-Type        Score          blaOXA-Genes")
        print(underscore + "__________________________________________________________")
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


def xspecT(BF, BF_1_1, files, paths, file_format, read_amount, metagenome):
    """performs a BF-lookup for a set of genomes for testing purpose"""
    print("Starting taxonomic assignment on species-level...")
    predictions = []
    scores = []
    counterx = 0
    for i in range(len(files)):
        if i == int(len(files)/6) or i == int(len(files)/3) or i == int(len(files)/2) or i == int(len(files)/1.5) or i == int(len(files)/1.2):
            print("...")
        BF.number_of_kmeres = 0
        BF.hits_per_filter = [0] * BF.clonetypes
        BF_1_1.number_of_kmeres = 0
        BF_1_1.hits_per_filter = [0]
        if file_format == "fasta" or file_format == "fna":
            for sequence in SeqIO.parse(paths[i], "fasta"):
                for j in range(0, len(sequence.seq) - BF.k, 500):
                    hits_per_filter_copy = BF_1_1.hits_per_filter[:]
                    BF_1_1.lookup(str(sequence.seq[j: j + BF.k]))
                    if hits_per_filter_copy != BF_1_1.hits_per_filter:
                        BF.number_of_kmeres += 1
                        BF.lookup(str(sequence.seq[j: j + BF.k]))
        elif file_format == "fastq" or file_format == "fq":
            counter = 0
            reads = []
            for sequence in SeqIO.parse(paths[i], "fastq"):
                #reverse_sequence = sequence.seq.reverse_complement()
                read_kmers = []
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
                if "N" not in str(sequence.seq):
                    BF_1_1.lookup(k1)
                    BF_1_1.lookup(k2)
                    BF_1_1.lookup(k3)
                    BF_1_1.lookup(k4)
                    BF_1_1.lookup(k5)
                # needs at least 2 of 3 hits to continue with read
                if (sum(BF_1_1.hits_per_filter) - hit_sum) > 3:
                    for j in range(len(str(sequence.seq)) - BF_1_1.k):
                        if "N" not in str(sequence.seq[j: j + BF_1_1.k]):
                            read_kmers.append(str(sequence.seq[j: j + BF_1_1.k]))
                    #read_kmers.append(str(reverse_sequence[j: j + BF_1_1.k]))
                    reads.append(read_kmers)
                    BF_1_1.hits_per_filter = hits_per_filter_copy
                else:
                    # resetting hit counter
                    BF_1_1.hits_per_filter = hits_per_filter_copy
            reads_filtered = set()
        #    reads_filtered = []
            threshold_dic = {}
            for i in range(len(reads)):
                read_kmers_filtered = []
                threshold = 0
                #hits_per_filter_copy = BF_1.hits_per_filter[:]
                #read_len = len(reads[i])
                #BF_1.lookup(reads[i][0])
                #BF_1.lookup(reads[i][int(read_len/2)])
                ##BF_1.lookup(reads[i][-1])
                #if (BF_1.hits_per_filter[0] - hits_per_filter_copy[0]) < 2:
                #    continue
                temp = []
                for j in range(len(reads[i])):
                    BF_1_1.number_of_kmeres += 1
                    hits_per_filter_copy = BF_1_1.hits_per_filter[:]
                    BF_1_1.lookup(reads[i][j])
                    if hits_per_filter_copy != BF_1_1.hits_per_filter:
                        threshold += 1
                        temp.append(reads[i][j])
                count = threshold_dic.get(threshold, 0)
                threshold_dic[threshold] = count + 1
                if threshold >= 80:
                    reads_filtered.update(temp)
                        #reads_filtered.add(reads[i][j])
                        #reads_filtered.append(reads[i][j])
                #reads_filtered.append(list(read_kmers_filtered_unique))
            x, y = [], []
            threshold_dic = dict(sorted(threshold_dic.items()))
            #for key, value in threshold_dic.items():
                #x.append(key)
                #y.append(value)
            #plt.plot(x, y)
            #plt.ylabel('frequency')
            #plt.xlabel('kmer-hits')
            #plt.title('Read-Hit Coverage')
            #plt.show()
            #plt.clf()
            print(dict(sorted(threshold_dic.items())))
            for i in reads_filtered:
                hits_per_filter_copy = BF.hits_per_filter[:]
                BF.lookup(i)
                if ((sum(BF.hits_per_filter) - sum(hits_per_filter_copy)) <= 5) and ((sum(BF.hits_per_filter) - sum(hits_per_filter_copy)) != 0):
                    BF.number_of_kmeres += 1
                elif (sum(BF.hits_per_filter) - sum(hits_per_filter_copy)) > 5:
                    BF.hits_per_filter = hits_per_filter_copy[:]
            #for i in range(len(reads_filtered)):
                #hits_per_filter_copy = BF.hits_per_filter[:]
                ##read_len = len(reads_filtered[i])
                #BF.lookup(reads[i][0])
                #BF.lookup(reads[i][int(read_len/2)])
                #BF.lookup(reads[i][-1])
                #for j in range(len(reads_filtered[i])):
                #    hits_per_filter_copy = BF.hits_per_filter[:]
                    #BF.lookup(reads_filtered[i][j])
                    #if (sum(BF.hits_per_filter) - sum(hits_per_filter_copy)) == 0:
                    #    continue
                    #if ((sum(BF.hits_per_filter) - sum(hits_per_filter_copy)) <= 5) and ((sum(BF.hits_per_filter) - sum(hits_per_filter_copy)) != 0):
                    #    BF.number_of_kmeres += 1
                        #if BF.hits_per_filter[50] != hits_per_filter_copy[50]:
                        #    print(reads_filtered[i][j])
                    #elif (sum(BF.hits_per_filter) - sum(hits_per_filter_copy)) > 5:
                    #    BF.hits_per_filter = hits_per_filter_copy[:]
            print("Kmers searched in Acinetobacter-Filter: ", BF_1_1.number_of_kmeres)
            print("Kmers found in Acinetobacter-Filter: ", BF_1_1.hits_per_filter[0])
            print("Kmers discarded: ", (BF_1_1.number_of_kmeres - BF_1_1.hits_per_filter[0]))
            print("Kmers used for species assignment: ", BF.number_of_kmeres)
        #    for sequence in SeqIO.parse(paths[i], "fastq"):
        #        if counter < read_amount:
        #            counter += 1
        #            for j in range(0, len(sequence.seq) - BF.k+1, 10):
        #                BF.number_of_kmeres += 1
        #                BF.lookup(str(sequence.seq[j: j + BF.k]))
        #        else:
        #            break
        #    counter = 0
    #        for sequence in SeqIO.parse(paths[i], "fastq"):
    #            reverse_sequence = sequence.seq.reverse_complement()
    #            if counter < read_amount:
    #                counter += 1
    #                for j in range(0, len(reverse_sequence) - BF.k+1, 10):
    #                    BF.number_of_kmeres += 1
    #                    BF.lookup(str(reverse_sequence[j: j + BF.k]))
    #            else:
    #                break

        score = BF.get_score()
        print("Scores: ", score)
        print("Hits: ", BF.hits_per_filter)
        names = []
        with open(r'filter/FilterSpecies.txt', 'rb') as fp:
            names = pickle.load(fp)
        score_edit = [str(x) for x in score]
        score_edit = ",".join(score_edit)
        # making prediction
        if not metagenome:
            prediction = Classifier.classify(r'Training_data/Training_data_spec.csv', score, True)
        else:
            index_result = max(range(len(score)), key=score.__getitem__)
            prediction = names[index_result]
        if max(score) < 0.3:
            predictions.append("unknown")
            scores.append(str(max(score)))
        else:
            predictions.append("A. " + prediction)
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
                #Originally 10
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
    paths_oxa = sorted(os.listdir(r"filter/OXAs/"))
    oxas = []
    scores_oxa = []
    for i in paths_oxa:
        oxas.append(i[:-4])

    for i in range(len(files)):
        oxa_dic = {}
        if i == int(len(files)/6) or i == int(len(files)/3) or i == int(len(files)/2) or i == int(len(files)/1.5) or i == int(len(files)/1.2):
            print("...")
        # Checking file type
        # if the file is fasta -> concat lines
        reads = []
        BF_3.number_of_kmeres = 0
        BF_3.hits_per_filter = [0] * BF_3.clonetypes
        BF_3.table = OXATable()
        BF_3.table.read_dic(r'filter/OXAs_dict/oxa_dict.txt')
        if file_format == "fasta" or file_format == "fna":
            for sequence in SeqIO.parse(paths[i], "fasta"):
                reads.append(str(sequence.seq))
            BF_3.lookup_oxa(reads, ".fna")
        elif file_format == "fastq" or file_format == "fq":
            counter = 0
            for sequence in SeqIO.parse(paths[i], "fastq"):
                if counter < read_amount:
                    counter += 1
                    reads.append(str(sequence.seq))
                else:
                    break
            BF_3.lookup_oxa(reads, ".fq")
            #print("Reads used: ", counter)
        score_oxa = BF_3.get_oxa_score()
        for i in range(len(oxas)):
            oxa_dic[oxas[i]] = score_oxa[i]
        for i in range(len(oxa_dic)):
            if oxa_dic[oxas[i]] < 0.3:
                del oxa_dic[oxas[i]]
        if len(oxa_dic) == 0:
            oxa_dic = "None"
        scores_oxa.append(oxa_dic)
    end = time.time()
    needed = round(end - start, 2)
    print("Time needed: ", needed)
    print("Screening for blaOXA-genes done...")
    return scores_oxa


def main():
    print("")
    print("XspecT performs a taxonomic assignment on the species-level for bacteria of the genus Acinetobacter.")
    print("ClAssT performs a strain-typing one sub-type-level for A. baumannii.")
    print("You can also screen your file for blaOXA-genes.")
    print("XspecT/ClAssT needs a file path to all files where an assignment shall be performed.")
    print("")
    print("Run XspecT: (y/n)?")
    XspecT = input()
    if XspecT == "y":
        XspecT = True
    elif XspecT == "n":
        XspecT = False
    else:
        print("Error: Wrong Input, use y/n (y=yes, n=no)!")
        quit()
    print("Run ClAssT: (y/n)?")
    ClAssT = input()
    if ClAssT == "y":
        ClAssT = True
    elif ClAssT == "n":
        ClAssT = False
    else:
        print("Error: Wrong Input, use y/n (y=yes, n=no)!")
        quit()
    print("Screen for blaOXA-Genes: (y/n)?")
    oxa = input()
    if oxa == "y":
        oxa = True
    elif oxa == "n":
        oxa = False
    else:
        print("Error: Wrong Input, use y/n (y=yes, n=no)!")
        quit()
    print("Metagenome-Mode?: (y/n)?")
    metagenome = input()
    if metagenome == "y":
        metagenome = True
    elif metagenome == "n":
        metagenome = False
    else:
        print("Error: Wrong Input, use y/n (y=yes, n=no)!")
        quit()
    print("Enter file-format: (fasta/fastq)")
    file_format = input()
    if file_format != "fasta" and file_format != "fastq":
        print("Error: Invalid file-format, use fasta/fastq!")
        quit()
    if file_format == "fastq":
        print("How many reads should be used: ")
        read_amount = int(input())
    else:
        read_amount = 543789
    if read_amount < 0:
        print("Error: Invalid read amount")
    print("Enter file-path:")
    file_path = input()
    print("Save results as CSV-table: (y/n)?")
    csv_table = input()
    print("")
    if csv_table == "y":
        csv_table = True
    elif csv_table == "n":
        csv_table = False
    else:
        print("Error: Wrong Input, use y/n (y=yes, n=no)!")
        quit()
    if XspecT == False and ClAssT == False and oxa == False:
        print("No tool selected, closing application...")
        quit()
    xspecT_mini(file_path, XspecT, ClAssT, oxa, file_format, read_amount, csv_table, metagenome)


if __name__ == '__main__':
    main()
