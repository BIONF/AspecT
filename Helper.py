import BF_v2
import search_filter
import pickle
from bitarray import bitarray
import os
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from Bio import SeqIO, SeqRecord, Seq
import csv
import Classifier
#from Bio.Seq import Seq

# Help-Script to generate new default Data, manually train new Bloomfilter, generate svm-data and test the tools

def write_file():
    # https://stackoverflow.com/questions/899103/writing-a-list-to-a-file-with-python/899176
    itemlist = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5','IC6','IC7','IC8']
    with open(r'filter/FilterClonetypes.txt', 'wb') as fp:
        pickle.dump(itemlist, fp)


def write_file2():
    # https://stackoverflow.com/questions/899103/writing-a-list-to-a-file-with-python/899176
    itemlist = ['bla OXA-23', 'bla OXA-24', 'bla OXA-51', 'bla OXA-58', 'bla OXA-134',
                'bla OXA-143','bla OXA-211','bla OXA-213','bla OXA-214','bla OXA-229','bla OXA-286']
    with open(r'filter/OXAs/FilterOXA.txt', 'wb') as fp:
        pickle.dump(itemlist, fp)


def write_file3():
    # https://stackoverflow.com/questions/899103/writing-a-list-to-a-file-with-python/899176
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
                "vivianii", "wanghuae", "wuhouensis"]

    with open(r'filter/FilterSpecies.txt', 'wb') as fp:
        pickle.dump(itemlist, fp)


def train():
    files = os.listdir(r'I:\OXA-Gene')
    for i in range(len(files) -1, -1, -1):
        if 'fasta' not in files[i]:
            del files[i]

    for i in range(len(files)):

        BF = BF_v2.AbaumanniiBloomfilter(80000)
        BF.set_arraysize(80000)
        BF.set_clonetypes(1)
        BF.set_hashes(7)
        BF.set_k(20)
        path = r'I:/OXA-Gene/' + files[i]
        name = files[i][:-6] + '.txt'
        print(name)
        result = r'C:/Users/SG/Desktop/' + name
        BF.train_sequence(path, 0)
        BF.save_clonetypes(result)
        BF.cleanup()


#changed directory
def train_Core():
    """trains (concatenated-)genomes into BF and saves them"""
    #files = os.listdir(r"filter\species_totrain")
    files = os.listdir(r'F:\project\genomes\totrain')
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
        #path = r"filter/species_totrain/" + files[i]
        path = r'F:/project/genomes/totrain/' + files[i]
        name = files[i].split('.')[-2] + '.txt'
        print(name)
        #result = r"filter/species" + name
        result = r'F:/project/results/' + name
        BF.train_sequence(path, 0)
        BF.save_clonetypes(result)
        BF.cleanup()


def opene():
    with open(r'C:\Users\SG\Desktop\a.baumannii Filter\FilterOXA.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    print(clonetypes)


def openspec():
    with open(r'C:\Users\Dominik\Uni\SoSe21\PraktikumBA-Arbeit\ClAssT-Acinetobacter-baumannii-Clone-type-Assignment-Tool-master\ClAssT-Acinetobacter-baumannii-Clone-type-Assignment-Tool-master\filter\FilterSpecies.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    for i in clonetypes:
        print(i)


def pw():
    from flask_bcrypt import Bcrypt
    bcrypt = Bcrypt()
    print(bcrypt.generate_password_hash('user'))
    print(bcrypt.generate_password_hash('pwd'))


def Test():
    temp = bitarray(0)
    with open(r'filter\OXA51_IC1.txt', 'rb') as fh:
        temp.fromfile(fh)
        print(len(temp))


def Test_Core_for_OXA():
    with open(r'filter/FilterClonetypes.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)

    BF = BF_v2.AbaumanniiBloomfilter(22000000)
    BF.set_arraysize(22000000)
    BF.set_hashes(7)
    BF.set_k(20)
    # User Options
    BF.set_reads(1000)

    paths = [r'filter/CoreIC1.txt',
             r'filter/CoreIC2.txt',
             r'filter/CoreIC3.txt',
             r'filter/CoreIC4.txt',
             r'filter/CoreIC5.txt',
             r'filter/CoreIC6.txt',
             r'filter/CoreIC7.txt',
             r'filter/CoreIC8.txt']

    BF.read_clonetypes(paths, clonetypes)

    Oxa_paths = [r'H:\bla-51-like\IC1\OXA69.fasta',
                 r'H:\bla-51-like\IC1\OXA92.fasta',
                 r'H:\bla-51-like\IC1\OXA107.fasta',
                 r'H:\bla-51-like\IC1\OXA110.fasta',
                 r'H:\bla-51-like\IC2\OXA66.fasta',
                 r'H:\bla-51-like\IC2\OXA82.fasta',
                 r'H:\bla-51-like\IC2\OXA172.fasta',
                 r'H:\bla-51-like\IC2\OXA201.fasta',
                 r'H:\bla-51-like\IC2\OXA202.fasta',
                 r'H:\bla-51-like\IC3\OXA71.fasta',
                 r'H:\bla-51-like\IC3\OXA113.fasta',
                 r'H:\bla-51-like\IC4\OXA51.fasta',
                 r'H:\bla-51-like\IC4\OXA219.fasta',
                 r'H:\bla-51-like\IC5\OXA65.fasta',
                 r'H:\bla-51-like\IC6\OXA90.fasta',
                 r'H:\bla-51-like\IC6\OXA200.fasta',
                 r'H:\bla-51-like\IC7\OXA64.fasta',
                 r'H:\bla-51-like\IC8\OXA68.fasta',
                 r'H:\bla-51-like\IC8\OXA128.fasta']

    for path in Oxa_paths:
        BF.lookup_sequence(path)
        score = BF.get_score()
        print(score)

def csv_helper():
    files = os.listdir(r'F:\project\test-set')
    for i in range(len(files) -1, -1, -1):
        if 'fna' in files[i] or 'fasta' in files[i]:
            continue
        else:
            del files[i]
    with open(r'F:/project/csv/help.csv', 'w') as file:
        writer = csv.writer(file)
        writer.writerows(files)

def distinct_kmer():
    """ creates a fasta-file with the distinct kmers of every species """
    files = os.listdir(r'F:\project\genomes\coverage')
    for i in range(len(files) -1, -1, -1):
        if 'fna' in files[i] or 'fasta' in files[i]:
            continue
        else:
            del files[i]
    paths = files[:]
    for i in range(len(files)):
        paths[i] = r'F:/project/genomes/coverage/' + paths[i]
    counter = 0
    for j in [41, 51]:
        for i in range(len(files)):
            counter += 1
            print(files[i])
            records = []
            kmers = []
            for sequence in SeqIO.parse(paths[i], "fasta"):
                for i in range(len(sequence.seq) - j + 1):
                    kmers.append(str(sequence.seq[i: i + j]))
            #print(len(kmers))
            distinct_kmer = []
            distinct_kmer = list(dict.fromkeys(kmers))
            print(len(distinct_kmer))
            for kmer in distinct_kmer:
                records.append(SeqRecord.SeqRecord(Seq.Seq(kmer)))
            with open(r"F:/project/genomes/coverage/result/distinct_complete_" + str(j) + "kmer_" + str(counter) + ".fasta", "a") as output_handle:
                SeqIO.write(records, r"F:/project/genomes/coverage/result/distinct_complete_" + str(j) + "kmer_" + str(counter) + ".fasta", "fasta")


def coverage_plot(min_coverage=7, max_coverage=1000):
    """creates a coverage plot from a histo-file, used for kmer-specific research"""
    #http://voorloopnul.com/blog/kmer-analysis-with-python-part-1/
    files = os.listdir(r'F:\project\genomes\coverage\result\results\histo')
    paths = files[:]
    for i in range(len(files)):
        paths[i] = r'F:/project/genomes/coverage/result/results/histo/' + paths[i]
    for i in range(len(files)):
        ext = files[i].split('.')[-2]
        if "14" in ext or "16" in ext or "18" in ext:
            continue
        else:
            with open(paths[i]) as file:
                data = file.readlines()

            dataset = [entry.replace("\n", "") for entry in data]
            dataset = [entry.split(" ") for entry in dataset]
            #print(dataset)
            coverage = [int(entry[0]) for entry in dataset][min_coverage:max_coverage]
            frequency = [int(entry[1]) for entry in dataset][min_coverage:max_coverage]
            higher_frequency = max(frequency)

            plt.plot(coverage, frequency, label="%s"%(ext[-6:],))
            leg = plt.legend(loc='lower right', ncol=1, shadow=True, fancybox=True)
            leg.get_frame().set_alpha(0.5)
            plt.ylabel('Frequency')
            plt.xlabel('Coverage')
            #plt.title(ext)
            plt.ylim(0, higher_frequency)
    plt.title("Acinetobacter-Coverage comparison with different kmer-lengths")
    plt.show()



def count_kmer():
    """creates multiple coverage plots from a set of genomes"""
    # https://stackoverflow.com/questions/2600191/how-can-i-count-the-occurrences-of-a-list-item
    # https://stackoverflow.com/questions/20950650/how-to-sort-counter-by-value-python
    # https://pythonexamples.org/python-sort-list-of-tuples/
    files = os.listdir(r'F:\project\genomes\totrain')
    for i in range(len(files) -1, -1, -1):
        if 'fna' in files[i] or 'fasta' in files[i]:
            continue
        else:
            del files[i]
    paths = files[:]
    for i in range(len(files)):
        paths[i] = r'F:/project/genomes/totrain/' + paths[i]
    for i in range(len(files)):
        print(files[i])
        ext = files[i].split('.')[-2]
        kmers = []
        intermediat_kmers = []
        for sequence in SeqIO.parse(paths[i], "fasta"):
            for i in range(len(sequence.seq) - 20 + 1):
                kmers.append(str(sequence.seq[i: i + 20]))
        kmers_freq = Counter(kmers)
        kmers_freq_freq = []
        kmers_no_duplicates = list(dict.fromkeys(kmers))
        for i in range(len(kmers_freq)):
            kmers_freq_freq.append(kmers_freq[kmers_no_duplicates[i]])
        kmers_freq_freq = Counter(kmers_freq_freq)
        kmers_freq_sorted = kmers_freq_freq.most_common()
        kmers_freq_sorted.sort(key = lambda x: x[0])
        kmer_index = 0
        print(kmers_freq_freq)
        print(kmers_freq_sorted)
        for i in range(1, len(kmers_freq_sorted)):
            if kmers_freq_sorted[i][1] > kmers_freq_sorted[i + 1][1]:
                kmer_index += 1
            else:
                break
        for y in kmers:
            if kmers_freq[y] > kmer_index:
                intermediat_kmers.append(y)
        intermediat_kmers = list(dict.fromkeys(intermediat_kmers))
        x_achse = []
        y_achse = []
        for i in range(len(kmers_freq_sorted)):
            x_achse.append(kmers_freq_sorted[i][0])
            y_achse.append(kmers_freq_sorted[i][1])
        plt.plot(x_achse[kmer_index - 1:-1], y_achse[kmer_index - 1:-1])
        plt.ylabel('frequency')
        plt.xlabel('kmer-coverage')
        plt.title('Distribution of '+ ext + '-kmers in Acinetobacter species')
        plt.savefig(r'F:/project/kmere-distr/kmer-coverage/full/' + (ext + ".png"))
        plt.clf()

def test_genomes():
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
    best_match = ["GCF_000587995",
                    "GCF_000588015",
                    "GCF_000588095",
                    "GCF_000588255",
                    "GCF_000588335",
                    "GCF_000588395",
                    "GCF_000588415",
                    "GCF_000588515",
                    "GCF_000588535",
                    "GCF_000588595",
                    "GCF_000588715",
                    "GCF_000589115",
                    "GCF_000589135",
                    "GCF_000589155",
                    "GCF_000589175",
                    "GCF_000589215",
                    "GCF_000589235",
                    "GCF_000589255",
                    "GCF_000589335",
                    "GCF_000620525",
                    "GCF_000681815",
                    "GCF_000682155",
                    "GCF_000682355",
                    "GCF_000682615",
                    "GCF_000705635",
                    "GCF_000788125",
                    "GCF_000800865",
                    "GCF_001005515",
                    "GCF_001005525",
                    "GCF_001054715",
                    "GCF_001077515",
                    "GCF_001276505",
                    "GCF_001278715",
                    "GCF_001423205",
                    "GCF_001425285",
                    "GCF_001471615",
                    "GCF_001541635",
                    "GCF_001592855",
                    "GCF_001605865",
                    "GCF_001720685",
                    "GCF_001729365",
                    "GCF_001866125",
                    "GCF_001878775",
                    "GCF_002018895",
                    "GCF_002018925",
                    "GCF_002251565",
                    "GCF_002412355",
                    "GCF_002797085",
                    "GCF_002797235",
                    "GCF_002797255",
                    "GCF_002899995",
                    "GCF_002918965",
                    "GCF_002919885",
                    "GCF_002934965",
                    "GCF_900110445",
                    "GCF_900110525"
                    ]
    print("Preparing BloomFilter...")
    BF = search_filter.pre_processing()
    print("Collecting Input-Data...")
    #files = os.listdir(r'F:\project\genomes\all_genomes_edit')
    files = os.listdir(r'F:\project\test-set')
    #files = os.listdir(r'F:\project\test-set\not_Acinetobacter')
    #files = os.listdir(r'F:\project\genomes\unclassified')
    for i in range(len(files) -1, -1, -1):
        if 'fna' in files[i] or 'fasta' in files[i]:
            continue
        else:
            del files[i]
    paths = files[:]
    for i in range(len(files)):
        #paths[i] = r'F:/project/genomes/all_genomes_edit/' + paths[i]
        paths[i] = r'F:/project/test-set/' + paths[i]
        #paths[i] = r'F:/project/test-set/not_Acinetobacter/' + paths[i]
        #paths[i] = r'F:/project/genomes/unclassified/' + paths[i]
    names = []
    print("Saving Species-Names...")
    for i in range(len(files)):
        with open(paths[i]) as file:
            head = file.readline()
            head = head.split()
            try:
                names.append(head[2])
            except:
                names.append("NameError")
    GCF_numbers = []
    print("Saving GCF_Numbers...")
    for i in range(len(files)):
        GCF = files[i].split('.')[0]
        GCF_numbers.append(GCF)
    print("Starting Taxonomic Assignment on Species-Level...")
    predictions = []
    scores = []
    #scores_plot
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
        scores.append(score_edit)
        prediction = Classifier.classify(r'Training_data/Training_data_spec+none_new.csv', score, False)
        predictions.append(prediction)
    #    if GCF_numbers[i] in best_match:
            #names[i] = prediction
        if names[i] in itemlist:
            test[itemlist.index(names[i])][itemlist.index(prediction)] += 1
        else:
            print(GCF_numbers[i])
            print("Falscher Name, aus Test ausgeschlossen")
        #scores_plot.append([name[i]] + [prediction] + score)
    #scores_plot = sorted(scores_plot, key = lambda h: h[0])
    #for i in range(len(files)):
        #for j in range(len(files)):
            #if scores_plot[i][0] = scores_plot[j][0]:
            #    ## TODO:
            #else:
            #    break
    for i in range(len(test)):
        summe = sum(test[i])
        for j in range(len(test[0])):
            try:
                test[i][j] = test[i][j] / summe
            except:
                continue
    print("Assignment Done...")
    print("Prepare Result-Data...")
    excel = []
    for i in range(len(files)):
        excel.append(GCF_numbers[i] + "," + scores[i] + "," + names[i] + "," + predictions[i])
    for i in range(len(excel)):
        excel[i] = [excel[i]]
    with open(r'F:/project/csv/Test-Set_20.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(excel)
    print("Finished!")


def test():
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
    x = [[2,1,1,1,],
        [1,1,5,1,],
        [1,1,1,1,],
        [1,1,1,1,],
        [1,1,1,3,]]
    test = [[0 for i in range(83)] for j in range(83)]
    plt.figure(figsize=(20,15))
    ax = sns.heatmap(test, xticklabels = itemlist, yticklabels = itemlist, cmap = "Blues", linewidth=1)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
    plt.show()


def main():
    #Test_Core_for_OXA()
    #write_file()
    #opene()
    #openspec()
    #write_file2()
    #pw()
    #train_Core()
    #write_file3()
    #write_file4()
    #histo()
    #count_kmer()
    #distinct_kmer()
    #coverage_plot()
    #csv_helper()
    test_genomes()
    #test()

if __name__ == '__main__':
    main()
