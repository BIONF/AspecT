from sklearn.svm import SVC
import csv
from copy import deepcopy
import pickle


def cut_csv(csv_file, lst, table=False):
    """Returns desired data from Training_data"""
    r = csv.reader(open(csv_file))
    m = list(r)
    selected = deepcopy(lst)

    header = m[0]
    m = m[1:]
    labels = header[1:-1]

    X_train = []
    y_train = []
    files = []

    if selected[8] and len(header[8:-1]) > 0:
        # Added Genomes selected
        del selected[8]
        selected = selected + ([True] * len(header[9:-1]))

    else:
        # Added Genomes not selected
        del selected[8]
        selected = selected + ([False] * len(header[9:-1]))

    # creating matrix
    for i in range(len(m)):
        X_train.append(m[i][1:-1])
        y_train.append(m[i][-1])
        files.append(m[i][0])

    # Deleting Cols
    for i in range(len(X_train)):
        for j in range(len(X_train[i]) - 1, -1, -1):
            if selected[j]:
                pass
            else:
                del X_train[i][j]


    # Deleting Rows
    valid = ['None']
    for i in range(len(selected)):
        if selected[i]:
            valid.append(labels[i])

    for i in range(len(X_train) - 1, -1, -1):
        if y_train[i] not in valid:
            del y_train[i]
            del X_train[i]
            del files[i]

    if table:
        # Inserting Infos for Table
        for i in range(len(X_train)):
            X_train[i].insert(0, files[i])
            X_train[i].append(y_train[i])

        for i in range(len(header) - 1, -1, -1):
            if header[i] not in valid:
                del header[i]

        header.insert(0, 'File')
        header.append('Label')

        X_train.insert(0, header)

    else:
        pass

    return X_train, y_train


def cut_csv_spec(csv_file):
    """Returns svm Training_data"""
    # read the training-data
    r = csv.reader(open(csv_file))
    m = list(r)
    header = m[0]
    m = m[1:]
    X_train = []
    y_train = []

    # creating matrix as input for the classifier
    for i in range(len(m)):
        X_train.append(m[i][1:-1])
        y_train.append(m[i][-1])

    return X_train, y_train


def classify(csv_file, result, lst):
    """ Classifys Result-vector and calculates needed vectors"""
    r = csv.reader(open(csv_file))
    m = list(r)
    # deciding which kernel-function will be used
    if m[0][1] == "IC1":
        mode = "ClAssT"
        X_train, y_train = cut_csv(csv_file, lst)
        svm = SVC(kernel='poly', C=1.0).fit(X_train, y_train)
    else:
        mode = "XspecT"
        X_train, y_train = cut_csv_spec(csv_file)
        svm = SVC(kernel='rbf', C=1.5).fit(X_train, y_train)
    # perform a prediction using the svm
    prediction = svm.predict([result])
    if mode == "XspecT":
        if max(result) < 0.3:
            prediction = ["sp.", 0]
    else:
        if max(result) < 0.3:
            prediction = ["None", 0]

    return prediction[0]


def IC3_classify(result_2):
    ic = 'International Clonetype 3 (ST32 or ST250)'
    m_3 = [['GCF_000278625.1', 1.0, ic],
           ['GCF_001674185.1', 0.86, ic],
           ['fictional', 0.85, 'NONE of the selected Clonetypes or Genomes'],
           ['fictional', 0.01, 'NONE of the selected Clonetypes or Genomes']]

    X = []
    y = []
    for i in range(len(m_3)):
        X.append(m_3[i][1])
        y.append(m_3[i][2])

    for i in range(len(X)):
        X[i] = [X[i]]
    svm_IC3 = SVC(kernel='poly', C=1).fit(X, y)

    return svm_IC3.predict([result_2]), result_2[0]

#https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
