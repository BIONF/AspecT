import csv
import os
import pickle
from linecache import getline
from pathlib import Path
from time import sleep, asctime, localtime

from Bio import SeqIO
from loguru import logger

import BF_v2
from train_filter.ncbi_api import download_assemblies
from train_filter import html_scrap


def get_current_time():
    """Returns the current time in the form hh:mm:ss."""
    return asctime(localtime()).split()[3]


def select_assemblies(accessions):
    """Selects up to 4 assemblies, preferably assemblies that were not used for training the filters.

    :param accessions: All selected assembly accessions for every species.
    :type accessions: dict
    :return: Dict with species name as key and selected accessions as value.
    """
    all_accessions = dict()

    for sci_name, current_accessions in accessions.items():
        selected_accessions = list()
        # Select 4 assemblies beginning from the last one.
        for i in range(len(current_accessions) - 1, -1, -1):
            selected_accessions.append(current_accessions[i])
            if len(selected_accessions) == 4:
                break

        all_accessions[sci_name] = selected_accessions

    return all_accessions


def get_svm_assemblies(all_accessions, dir_name):
    """Download assemblies for svm creation.

    :param all_accessions: Contains lists with all previously selected assemblies for every species.
    :type all_accessions: dict
    :param dir_name: Name of the parent directory.
    :type dir_name: str
    """
    # Select accessions for download.
    selected_accessions = select_assemblies(all_accessions)

    # Download assemblies.
    for sci_name, accessions in selected_accessions.items():
        sleep(5)
        logger.info("Downloading {name}", name=sci_name)
        file_name = sci_name + ".zip"
        download_assemblies.download_assemblies(
            accessions=accessions,
            dir_name=dir_name,
            target_folder="training_data_zipped",
            zip_file_name=file_name
        )
    logger.info("Downloads finished")


def delete_non_fasta(files):
    """Delete all non fasta files from the list.

    :param files: List of file names.
    :type files: list[str]
    :return: List with only fasta files.
    """
    # All possible fasta file endings.
    fasta_endings = ["fasta", "fna", "fa", "ffn", "frn"]

    # Iterate through file list backwards and delete all non fasta files.
    for i in range(len(files) - 1, -1, -1):
        file = files[i].split(".")
        if file[-1] in fasta_endings:
            continue
        else:
            del files[i]

    return files


def get_accessions(files):
    """Extracts the accession from the file names.

    :param files: List of file names.
    :type files: list[str]
    :return: List of all accessions.
    """
    accessions = list()
    for i in range(len(files)):
        accessions.append(files[i].split("_"))
        accessions[i] = accessions[i][0] + "_" + accessions[i][1]

    return accessions


def get_paths(dir_name, files):
    """Make a list with the paths to the files.

    :param dir_name: Name of the parent directory.
    :type dir_name: str
    :param files: List of file names.
    :type files: list[str]
    :return: A list with all file paths.
    """
    path = Path(__file__).parent.parent / "genus_metadata" / dir_name / "training_data"
    paths = list()
    for file in files:
        paths.append(str(path / file))

    return paths


def get_species_names(file_paths):
    """Extracts the species names.

    :param file_paths: List with the file paths.
    :type file_paths: list[str]
    :return: List with all species names.
    """
    names = list()
    for path in file_paths:
        header = getline(path, 1)
        name = header.replace("\n", "").replace(">", "")
        if not name.isdigit():
            logger.error("The header of file: {path} does not contain a correct ID: {name}. The ID needs to be "
                         "just numbers")
            logger.error("Aborting")
            exit()
        names.append(name)
    return names


def init_bf(genus, array_size, hashes=7, k=21):
    """Initializes bloomfilter.

    :param genus: Name of the genus.
    :type genus: str
    :param array_size: Size of the bloomfilter.
    :type array_size: int
    :param hashes: The number of hash functions the bf uses.
    :type hashes: int
    :param k: Length of k-mers.
    :type k: int
    :return: The bloomfilter object.
    """
    path = Path(__file__).parent.parent / "filter"

    # Initialize bloomfilter for genus.
    BF = BF_v2.AbaumanniiBloomfilter(array_size)
    BF.set_arraysize(array_size)
    BF.set_hashes(hashes)
    BF.set_k(k)

    # Get all species names.
    names_path = path / "species_names" / ("Filter" + genus + ".txt")
    with open(names_path, "rb") as f:
        clonetypes = pickle.load(f)

    # Get bloomfilter paths.
    bf_path = path / genus
    paths = sorted(os.listdir(bf_path))
    for i in range(len(paths)):
        paths[i] = str(bf_path / str(paths[i]))
    # Setup bloomfilters.
    BF.read_clonetypes(paths, clonetypes)

    return BF


def perform_lookup(bloomfilter, files, file_paths, accessions, names, spacing):
    """Performs a lookup on a bloomfilter object and gives the scores as a list.

    :param bloomfilter: The bloomfilter object on which the lookup is performed.
    :param files: List of file names.
    :type files: list[str]
    :param file_paths: List with the file paths.
    :type file_paths: list[str]
    :param accessions: List of all accessions.
    :type accessions: list[str]
    :param names: List with all species names.
    :type names: list[str]
    :return: List with all scores of the lookup.
    """
    scores = list()
    BF = bloomfilter

    # Lookup.
    for i in range(len(files)):
        BF.number_of_kmeres = 0
        BF.hits_per_filter = [0] * BF.clonetypes

        for sequence in SeqIO.parse(file_paths[i], "fasta"):
            # Dominik: changed sample size to var
            for j in range(0, len(sequence.seq) - BF.k, spacing):
                BF.number_of_kmeres += 1
                BF.lookup(str(sequence.seq[j: j + BF.k]))

        score = BF.get_score()
        score = [str(x) for x in score]
        score = ",".join(score)
        scores.append(accessions[i] + "," + score + "," + names[i])

    return scores


# https://stackoverflow.com/questions/21431052/sort-list-of-strings-by-a-part-of-the-string
def sort_list(scores, names):
    """Sorts the scores list by species name.

    :param scores: The scores gathered by a lookup of a bloomfilter.
    :type scores: list
    :param names: List with all species names.
    :type names: list[str]
    :return: The sorted scores list.
    """
    scores.sort(key=lambda x: x.split(",")[-1][:2])
    names = [x for x in names if x != "none"]
    names = list(dict.fromkeys(names))
    scores.insert(0, sorted(names))
    scores[0] = ["File"] + scores[0] + ["Label"]

    for i in range(1, len(scores)):
        line = scores[i].split(",")
        scores[i] = line

    return scores


def save_csv(genus, scores):
    """Saves the scores as csv file.

    :param genus: Name of the genus.
    :type genus: str
    :param scores: The scores gathered by a lookup of a bloomfilter.
    :type scores: list
    """
    path = Path(__file__).parent.parent / "Training_data" / (genus + "_Training_data_spec.csv")
    with open(path, "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(scores)


# Dominik: added spacing
def new_helper(spacing, genus, dir_name, array_size, k=21):
    """Create support vector machine for bloomfilters of a genus.

    :param spacing:
    :param genus: Name of the genus.
    :type genus: str
    :param dir_name: Name of the parent directory.
    :type dir_name: str
    :param array_size: Size for the byte array which is the bloomfilter.
    :type array_size: int
    :param k: Length of the k-mers.
    :type k: int
    """
    # Get all files.
    path = Path(__file__).parent.parent / "genus_metadata" / dir_name / "training_data"
    files = os.listdir(path)

    # Delete all non fasta files.
    files = delete_non_fasta(files)

    # Get accessions from file names.
    accessions = get_accessions(files)

    # Get all complete file paths.
    file_paths = get_paths(dir_name, files)

    # Get all species names from the header in the fasta files.
    names = get_species_names(file_paths)

    # Initialize bloomfilter.
    bf = init_bf(genus, array_size)

    # Perform lookup on bloomfilter.
    # Dominik: added spacing
    scores = perform_lookup(bf, files, file_paths, accessions, names, spacing)

    # Sort score list by species names.
    scores = sort_list(scores, names)

    # Save results in csv file.
    save_csv(genus, scores)


def main():
    taxons = ['54736', '28901']
    used_accessions = ['GCF_000439255.1', 'GCF_006051015.1', 'GCF_006113225.1', 'GCF_016653555.1', 'GCF_000006945.2',
                       'GCF_000783815.2', 'GCF_001558355.2', 'GCF_001647755.1']
    array_size = 67000000
    dir_name = "Salmonella_15_9_2023_10-59-31"
    genus = "Salmonella"
    ani = html_scrap.TaxonomyCheck()
    ani_gcf = ani.ani_gcf
    print(f"Training new svm")
    new_helper(genus, dir_name, array_size, k=21)


if __name__ == "__main__":
    main()
