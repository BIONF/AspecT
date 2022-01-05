import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVC
import csv
from sklearn.metrics import confusion_matrix
from sklearn.metrics import plot_confusion_matrix
import random
from sklearn.model_selection import StratifiedKFold

def train_and_plot(kernel, X_train, y_train, X_test, y_test, title):
    svm = SVC(kernel=kernel, C=1).fit(X_train, y_train)
    svm_predictions = svm.predict(X_test)

    # model accuracy for X_test
    accuracy = str(round(svm.score(X_test, y_test), 4))
    title = title + accuracy
    # creating a confusion matrix
    cm = confusion_matrix(y_test, svm_predictions)

    # Plot non-normalized confusion matrix
    class_names = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'IC6', 'IC7', 'IC8', 'None']
    disp = plot_confusion_matrix(svm, X_test, y_test,
                                 display_labels=class_names,
                                 cmap=plt.cm.Blues,
                                 normalize='true',
                                values_format='.3f')
    # Source:
    # https://stackoverflow.com/questions/57043260/how-change-the-color-of-boxes-in-confusion-matrix-using-sklearn

    disp.ax_.set_title(title)
    plt.show()


def read_files(file):
    """Reads csv file and extracts necessary data for plotting"""
    r = csv.reader(open(file))
    lines = list(r)
    return lines

def get_ICs():
    """reads in csv with labels and genomes, assigns them to dictionary"""
    # Only taking genomes that have been assigned by MLST
    r = csv.reader(open("genomes_and_labels.csv"))
    lines = list(r)
    ICs = {
        "IC1": [],
        "IC2": [],
        "IC3": [],
        "IC4": [],
        "IC5": [],
        "IC6": [],
        "IC7": [],
        "IC8": [],
        "None": []
    }

    for i in range(len(lines)):
        labels = [lines[i][2], lines[i][4]]
        if labels[0] == labels[1]:
            # Oxford and Pasteur MLST assigned same label
            ICs[labels[0]].append(lines[i][0])
        else:
            # Using the protocol, that classified an IC
            # Could be either Oxford or Pasteur
            if labels[0] == 'None':
                del labels[0]
            ICs[labels[0]].append(lines[i][0])

    return ICs


def get_scores(file, dic, keys):
    r = csv.reader(open(file))
    lines = list(r)

    for key in keys:
        for i in lines:
            if i[0] in dic[key]:
                pos = dic[key].index(i[0])
                # Converting data to float
                dic[key][pos] = [float(x) for x in i[1:]]
    return dic


def get_random_training_vectors(dic, keys, n):
    # Chooses n random vectors per class, returns them
    # as np.arrays and separates training/test-data
    X_train = []
    y_train = []
    X_test = []
    y_test = []

    # Choosing random training vectors and removing them from test-data
    for key in keys:
        pos = random.sample(range(0, len(dic[key])), n)
        pos.sort()
        for i in range(len(pos)-1, -1, -1):
            X_train.append(dic[key][pos[i]])
            del dic[key][pos[i]]
            y_train.append(key)
        # Creating test-data
        X_test = X_test + dic[key]
        y_test = y_test + [key] * len(dic[key])

    return np.array(X_train), np.array(y_train), np.array(X_test), np.array(y_test)


def get_read_data(file, dic, keys):
    """ Assigns reads to labels"""
    r = csv.reader(open(file))
    lines = list(r)
    vecs_forwards = []
    labels_forwards = []
    vecs_reverse = []
    labels_reverse = []
    for key in keys:
        for i in dic[key]:
            for j in lines:
                if i in j[0]:
                    if '_2.fq' in j[0] or '_R2_' in j[0]:
                        vecs_reverse.append(j[2:])
                        labels_reverse.append(key)
                    else:
                        vecs_forwards.append(j[2:])
                        labels_forwards.append(key)

    return np.array(vecs_forwards), np.array(labels_forwards), np.array(vecs_reverse), np.array(labels_reverse)


def print_stats(acc):
    print("lowest Accuracy: ", min(acc))
    print("highest Accuracy: ", max(acc))
    print("average Accuracy: ", round((sum(acc) / len(acc)), 4))
    print(" ")
    print(" ")


def get_cm_and_accuracy(X_train, y_train, X_test, y_test, kernel, matrix=False):
    svm = SVC(kernel=kernel, C=1).fit(X_train, y_train)
    svm_predictions = svm.predict(X_test)

    # model accuracy for X_test
    accuracy = str(round(svm.score(X_test, y_test), 4))

    # creating a confusion matrix
    if matrix:
        print('Accuracy: ', accuracy)
        cm = confusion_matrix(y_test, svm_predictions)
        print(cm)
        print(' ')

    return float(accuracy)


def create_plots():
    ICs_assembly = get_ICs()
    ICs_core = get_ICs()

    # Changing NCBI accession to Score-Vector
    keys = ["IC1", "IC2", "IC3", "IC4", "IC5", "IC6", "IC7", "IC8", "None"]

    # Full Run: all k-mers are tested
    ICs_assembly_full = get_scores("Assemblys as reference, all k-mers.csv", ICs_assembly, keys)
    ICs_core_full = get_scores("Core-Genome as reference, all k-mers.csv", ICs_core, keys)
    # Quick-run: only every 10th k-mer is used
    ICs_assembly_quick = get_scores("Assemblys as reference, quick run.csv", ICs_assembly, keys)
    ICs_core_quick = get_scores("Core-Genome as reference, quick run.csv", ICs_core, keys)

    assembly_full = []
    assembly_full_label = []
    assembly_quick = []
    assembly_quick_label = []
    core_full = []
    core_full_label = []
    core_quick = []
    core_quick_label = []

    for key in keys:
        # Assembly Reference - full
        assembly_full = assembly_full + ICs_assembly_full[key]
        assembly_full_label = assembly_full_label + [key] * len(ICs_assembly_full[key])
        # Assembly Reference - quick
        assembly_quick = assembly_quick + ICs_assembly_quick[key]
        assembly_quick_label = assembly_quick_label + [key] * len(ICs_assembly_quick[key])

        # Core Reference - full
        core_full = core_full + ICs_core_full[key]
        core_full_label = core_full_label + [key] * len(ICs_core_full[key])
        # Assembly Reference - quick
        core_quick = core_quick + ICs_core_quick[key]
        core_quick_label = core_quick_label + [key] * len(ICs_core_quick[key])

    assembly_full, assembly_full_label = np.array(assembly_full), np.array(assembly_full_label)
    core_full, core_full_label = np.array(core_full), np.array(core_full_label)
    skf = StratifiedKFold(n_splits=10)
    skf.get_n_splits(assembly_full, assembly_full_label)

    print('Please Note: Data from Quick-run (only every 10th k-mer is tested) is used for testing only, meaning that the'
          'trainingdata is taken from the k-folded Full-run, while the quick-run is used for testing only.')
    print(' ')

    # Trying all possible Combinations

    # Assembly reference
    print(' Assembled genomes as reference ')
    print('StratifiedKFold (n=10): Assembled Genomes as reference, all k-mers are tested, polynomial SVM-Kernel')
    acc = []
    for train_index, test_index in skf.split(assembly_full, assembly_full_label):
        X_train, X_test = assembly_full[train_index], assembly_full[test_index]
        y_train, y_test = assembly_full_label[train_index], assembly_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, X_test, y_test, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Assembled Genomes as reference, all k-mers are tested, linear SVM-Kernel')
    acc = []
    for train_index, test_index in skf.split(assembly_full, assembly_full_label):
        X_train, X_test = assembly_full[train_index], assembly_full[test_index]
        y_train, y_test = assembly_full_label[train_index], assembly_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, X_test, y_test, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Assembled Genomes as reference, only every 10th k-mer is tested, polynomial SVM-Kernel')
    acc = []
    for train_index, test_index in skf.split(assembly_full, assembly_full_label):
        X_train, X_test = assembly_full[train_index], assembly_full[test_index]
        y_train, y_test = assembly_full_label[train_index], assembly_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, assembly_quick, assembly_quick_label, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Assembled Genomes as reference, only every 10th k-mer is tested, linear SVM-Kernel')
    acc = []
    for train_index, test_index in skf.split(assembly_full, assembly_full_label):
        X_train, X_test = assembly_full[train_index], assembly_full[test_index]
        y_train, y_test = assembly_full_label[train_index], assembly_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, assembly_quick, assembly_quick_label, 'linear'))
    print_stats(acc)

    # Core Genome
    print(' Core genome as reference ')
    print('StratifiedKFold (n=10): Core Genome as reference, all k-mers are tested, polynomial SVM-Kernel')
    acc = []
    for train_index, test_index in skf.split(core_full, core_full_label):
        X_train, X_test = core_full[train_index], core_full[test_index]
        y_train, y_test = core_full_label[train_index], core_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, X_test, y_test, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Core Genome as reference, all k-mers are tested, linear SVM-Kernel')
    acc = []
    for train_index, test_index in skf.split(core_full, core_full_label):
        X_train, X_test = core_full[train_index], core_full[test_index]
        y_train, y_test = core_full_label[train_index], core_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, X_test, y_test, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Core Genome as reference, only every 10th k-mer is tested, polynomial SVM-Kernel')
    acc = []
    for train_index, test_index in skf.split(core_full, core_full_label):
        X_train, X_test = core_full[train_index], core_full[test_index]
        y_train, y_test = core_full_label[train_index], core_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, core_quick, core_quick_label, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Core Genome as reference, only every 10th k-mer is tested, linear SVM-Kernel')
    acc = []
    for train_index, test_index in skf.split(core_full, core_full_label):
        X_train, X_test = core_full[train_index], core_full[test_index]
        y_train, y_test = core_full_label[train_index], core_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, core_quick, core_quick_label, 'linear'))
    print_stats(acc)

    # Testing Sequence Reads
    ICs_assembly = get_ICs()
    ICs_core = get_ICs()
    forward, forward_label, reverse, reverse_label = get_read_data('Reads as Input, Assembly as Reference.csv',
                                                                   ICs_assembly, keys)

    print('StratifiedKFold (n=10): Assembled genomes as reference, unmodified forward reads')
    acc = []
    for train_index, test_index in skf.split(assembly_full, assembly_full_label):
        X_train, X_test = assembly_full[train_index], assembly_full[test_index]
        y_train, y_test = assembly_full_label[train_index], assembly_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, forward, forward_label, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Assembled genomes as reference, unmodified reverse reads')
    acc = []
    for train_index, test_index in skf.split(assembly_full, assembly_full_label):
        X_train, X_test = assembly_full[train_index], assembly_full[test_index]
        y_train, y_test = assembly_full_label[train_index], assembly_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, reverse, reverse_label, 'poly'))

    print_stats(acc)

    forward, forward_label, reverse, reverse_label = get_read_data('Reads as Input, Assembly as Reference, 15% of the nucleotides changed.csv',
                                                                   ICs_assembly, keys)
    print('StratifiedKFold (n=10): Assembled genomes as reference, forward reads with 15% randomly changed nucleotides')
    acc = []
    for train_index, test_index in skf.split(assembly_full, assembly_full_label):
        X_train, X_test = assembly_full[train_index], assembly_full[test_index]
        y_train, y_test = assembly_full_label[train_index], assembly_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, forward, forward_label, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Assembled genomes as reference, reverse reads with 15% randomly changed nucleotides')
    acc = []
    for train_index, test_index in skf.split(assembly_full, assembly_full_label):
        X_train, X_test = assembly_full[train_index], assembly_full[test_index]
        y_train, y_test = assembly_full_label[train_index], assembly_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, reverse, reverse_label, 'poly'))

    print_stats(acc)

    # Core Genome Reference
    forward, forward_label, reverse, reverse_label = get_read_data('Reads as Input, Assembly as Reference.csv',
                                                                   ICs_core, keys)
    print('StratifiedKFold (n=10): Core Genome as reference, unmodified forward reads')
    acc = []
    for train_index, test_index in skf.split(core_full, core_full_label):
        X_train, X_test = core_full[train_index], core_full[test_index]
        y_train, y_test = core_full_label[train_index], core_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, forward, forward_label, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Core Genome as reference, unmodified reverse reads')
    acc = []
    for train_index, test_index in skf.split(core_full, core_full_label):
        X_train, X_test = core_full[train_index], core_full[test_index]
        y_train, y_test = core_full_label[train_index], core_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, reverse, reverse_label, 'poly'))

    print_stats(acc)


    forward, forward_label, reverse, reverse_label = get_read_data(
        'Reads as Input, Core-Genome as Reference, 15% of the nucleotides changed.csv',
        ICs_assembly, keys)
    print('StratifiedKFold (n=10): Core Genome as reference, forward reads with 15% randomly changed nucleotides')
    acc = []
    for train_index, test_index in skf.split(core_full, core_full_label):
        X_train, X_test = core_full[train_index], core_full[test_index]
        y_train, y_test = core_full_label[train_index], core_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, forward, forward_label, 'poly'))

    print_stats(acc)

    print('StratifiedKFold (n=10): Core Genome as reference, reverse reads with 15% randomly changed nucleotides')
    acc = []
    for train_index, test_index in skf.split(core_full, core_full_label):
        X_train, X_test = core_full[train_index], core_full[test_index]
        y_train, y_test = core_full_label[train_index], core_full_label[test_index]

        acc.append(get_cm_and_accuracy(X_train, y_train, reverse, reverse_label, 'poly'))

    print_stats(acc)

create_plots()





