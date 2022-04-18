import search_filter
import os
from Bio import SeqIO, SeqRecord, Seq
import Classifier



def test_genomes(file_path):
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
    BF = search_filter.pre_processing()
    files = os.listdir(file_path)
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
            print("Unkown, perhabs an unsupported species")
    for i in range(len(test)):
        summe = sum(test[i])
        for j in range(len(test[0])):
            try:
                test[i][j] = test[i][j] / summe
            except:
                continue
    print("Assignment done...")
    print("Preparing results...")
    print("")
    header_filename = "Filename"
    spaces = []
    space = "           "
    underscore = "________"
    name_max = len(max(itemlist, key=len))
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
    for i in range(len(files)):
        excel.append(files[i] + spaces[i] + predictions[i] + "       " + scores[i])
    print(header_filename + "           Prediction               Score")
    print(underscore + "_________________________________________")
    for i in excel:
        print(i)
    print("")
    print("")

def main():
    print("AspecT performs a taxonomic assignment on the species-level for bacteria of the genus Acinetobacter.")
    print("AspecT needs a file path to all files where an assignment shall be performed.")
    print("Enter file-path:")
    file_path = input()
    test_genomes(file_path)


if __name__ == '__main__':
    main()
