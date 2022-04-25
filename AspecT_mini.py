import search_filter
import os
from Bio import SeqIO, SeqRecord, Seq
import Classifier
from OXA_Table import OXATable
import warnings

warnings.filterwarnings("ignore")

def xspecT_mini(file_path, XspecT, ClAssT, oxa):
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
    if XspecT:
        BF = search_filter.pre_processing()
    if ClAssT:
        BF_2 = search_filter.pre_processing_ClAssT()
    if oxa:
        BF_3 = search_filter.pre_processing_oxa()
    try:
        files = sorted(os.listdir(file_path))
    except FileNotFoundError:
        print("Error: Invalid filepath!")
        quit()
    for i in range(len(files) -1, -1, -1):
        if 'fna' in files[i] or 'fasta' in files[i]:
            continue
        else:
            del files[i]
    paths = files[:]
    file_path2 = file_path[:]
    for i in range(len(file_path2)):
        if file_path2[i] == "\\":
            list_temp = list(file_path2)
            list_temp[i] = '/'
            file_path2 = ''.join(list_temp)
    for i in range(len(files)):
        paths[i] = file_path2 + "/" + paths[i]
    names = []
    for i in range(len(files)):
        with open(paths[i]) as file:
            head = file.readline()
            head = head.split()
            try:
                names.append(head[2])
            except:
                names.append("NameError")
    GCF_numbers = []
    for i in range(len(files)):
        GCF = files[i].split('.')[0]
        GCF_numbers.append(GCF)
    if XspecT:
        predictions, scores = xspecT(BF, files, paths, names, GCF_numbers)
    if ClAssT:
        predictions_ClAssT, scores_ClAssT = clAssT(BF_2, files, paths)
    if oxa:
        scores_oxa = blaOXA(BF_3, files, paths)
    print("Preparing results...")
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
    if ClAssT:
        for i in range(len(predictions_ClAssT)):
            if predictions_ClAssT[i] != "none":
                predictions_ClAssT[i] += " "
    if XspecT and ClAssT:
        for i in range(len(scores_ClAssT)):
            if scores[i] == "1.0":
                scores[i] += " "

    if XspecT and ClAssT and oxa:
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i] + "           " + predictions_ClAssT[i] + "            " + scores_ClAssT[i] + "           " + str(scores_oxa[i]))
        print(header_filename + "           Species                  Score          Sub-Type        Score          blaOXA-Genes")
        print(underscore + "__________________________________________________________________________________________________")
        for i in excel:
            print(i)
        print("")
        print("")
    elif XspecT and not ClAssT and not oxa:
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i])
        print(header_filename + "           Species                  Score")
        print(underscore + "_________________________________________")
        for i in excel:
            print(i)
        print("")
        print("")
    elif ClAssT and not XspecT and not oxa:
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] + predictions_ClAssT[i] + "            " + scores_ClAssT[i])
        print(header_filename + "           Sub-Type        Score")
        print(underscore + "________________________________")
        for i in excel:
            print(i)
        print("")
        print("")
    elif oxa and not ClAssT and not XspecT:
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] + str(scores_oxa[i]))
        print(header_filename + "           blaOXA-Genes")
        print(underscore + "___________________________")
        for i in excel:
            print(i)
        print("")
        print("")
    elif XspecT and ClAssT and not oxa:
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i] + "           " + predictions_ClAssT[i] + "            " + scores_ClAssT[i])
        print(header_filename + "           Species                  Score          Sub-Type        Score")
        print(underscore + "________________________________________________________________________")
        for i in excel:
            print(i)
        print("")
        print("")
    elif XspecT and oxa and not ClAssT:
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i]  + "           " + str(scores_oxa[i]))
        print(header_filename + "           Species                  Score          blaOXA-Genes")
        print(underscore + "___________________________________________________________________")
        for i in excel:
            print(i)
        print("")
        print("")
    elif ClAssT and oxa and not XspecT:
        for i in range(len(files)):
            excel.append(files[i] + spaces[i] + predictions_ClAssT[i] + "            " + scores_ClAssT[i] + "           " + str(scores_oxa[i]))
        print(header_filename + "           Sub-Type        Score          blaOXA-Genes")
        print(underscore + "__________________________________________________________")
        for i in excel:
            print(i)
        print("")
        print("")


def xspecT(BF, files, paths, names, GCF_numbers):
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
    print("Starting taxonomic assignment on species-level...")
    predictions = []
    scores = []
    test = [[0 for i in range(83)] for j in range(83)]
    for i in range(len(files)):
        if i == int(len(files)/6) or i == int(len(files)/3) or i == int(len(files)/2) or i == int(len(files)/1.5) or i == int(len(files)/1.2):
            print("...")
        BF.number_of_kmeres = 0
        BF.hits_per_filter = [0] * BF.clonetypes
        for sequence in SeqIO.parse(paths[i], "fasta"):
            for j in range(0, len(sequence.seq) - BF.k, 500):
                BF.number_of_kmeres += 1
                BF.lookup(str(sequence.seq[j: j + BF.k]))
        score = BF.get_score()
        score_edit = [str(x) for x in score]
        score_edit = ",".join(score_edit)
        prediction = Classifier.classify(r'Training_data/Training_data_spec.csv', score, False)
        predictions.append("A. " + prediction)
        scores.append(str(max(score)))
        if names[i] in itemlist:
            test[itemlist.index(names[i])][itemlist.index(prediction)] += 1
        else:
            print(GCF_numbers[i])
            print("Unkown, perhabs a unsupported species")
    for i in range(len(test)):
        summe = sum(test[i])
        for j in range(len(test[0])):
            try:
                test[i][j] = test[i][j] / summe
            except:
                continue
    print("Taxonomic assignment done...")
    return predictions, scores


def clAssT(BF_2, files, paths):
    print("Starting strain-typing on sub-type-level...")
    predictions_ClAssT = []
    scores_ClAssT = []
    for i in range(len(files)):
        if i == int(len(files)/6) or i == int(len(files)/3) or i == int(len(files)/2) or i == int(len(files)/1.5) or i == int(len(files)/1.2):
            print("...")
        BF_2.number_of_kmeres = 0
        BF_2.hits_per_filter = [0] * BF_2.clonetypes
        for sequence in SeqIO.parse(paths[i], "fasta"):
            for j in range(0, len(sequence.seq) - BF_2.k, 10):
                BF_2.number_of_kmeres += 1
                BF_2.lookup(str(sequence.seq[j: j + BF_2.k]))
        score_ClAssT = BF_2.get_score()
        score_edit_ClAssT = [str(x) for x in score_ClAssT]
        score_edit_ClAssT = ",".join(score_edit_ClAssT)
        prediction_ClAssT = Classifier.classify(r'Training_data/Training_data_IC.csv', score_ClAssT, [True,True,True,True,True,True,True,True,False])
        predictions_ClAssT.append(prediction_ClAssT)
        scores_ClAssT.append(str(max(score_ClAssT)))

    print("Strain-typing on sub-type-level done...")
    return predictions_ClAssT, scores_ClAssT


def blaOXA(BF_3, files, paths):
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
        for sequence in SeqIO.parse(paths[i], "fasta"):
            reads.append(str(sequence.seq))
        BF_3.number_of_kmeres = 0
        BF_3.hits_per_filter = [0] * BF_3.clonetypes
        BF_3.table = OXATable()
        BF_3.table.read_dic(r'filter/OXAs_dict/oxa_dict.txt')
        BF_3.lookup_oxa(reads, ".fna")
        score_oxa = BF_3.get_oxa_score()
        for i in range(len(oxas)):
            oxa_dic[oxas[i]] = score_oxa[i]
        for i in range(len(oxa_dic)):
            if oxa_dic[oxas[i]] < 0.3:
                del oxa_dic[oxas[i]]
        if len(oxa_dic) == 0:
            oxa_dic = "None"
        scores_oxa.append(oxa_dic)
    print("Screening for blaOXA-genes done...")
    return scores_oxa


def main():
    print("")
    print("XspecT performs a taxonomic assignment on the species-level for bacteria of the genus Acinetobacter.")
    print("ClAssT performs a strain-typing one sub-type-level for A. baumannii.")
    print("You can also screen your file for blaOXA-genes.")
    print("XspecT/ClAssT needs a file path to all files where an assignment shall be performed.")
    print("")
    print("Run Xspect: (y/n)?")
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
    print("Enter file-path:")
    file_path = input()
    print("")
    if XspecT == False and ClAssT == False and oxa == False:
        print("No tool selected, closing application...")
        quit()
    xspecT_mini(file_path, XspecT, ClAssT, oxa)


if __name__ == '__main__':
    main()
