import subprocess as sp
from linecache import getline
from os import listdir, remove
from pathlib import Path
from time import perf_counter, localtime, asctime

from loguru import logger
import numpy as np

import XspecT_trainer


def get_seq_paths(dir_name: str):
    """Stores the sequence paths, with the species name as key, in a dictionary. The sequences are DNA-Assemblies which
    were concatenated.

    :param dir_name: Name of the directory.
    :return: Dictionary with species names and sequences.
    """
    dir_path = Path(__file__).parent.parent / 'genus_metadata' / dir_name / 'concatenate'
    sequence_dict = dict()
    files = listdir(dir_path)
    # Go through all files backwards to delete all non fasta files.
    for i in range(len(files) -1, -1, -1):
        curr_file = str(files[i])
        file_parts = curr_file.split(".")
        if file_parts[-1] != "fasta":
            del files[i]
        else:
            species_name = file_parts[0]
            # Save the species name with the path to its fasta file.
            sequence_dict[species_name] = str(dir_path / curr_file)

    return sequence_dict


def jellyfish_count(command: str):
    """A jellyfish command to count the k-mers of an fasta file using the linux bash.

    :param command: The jellyfish command with all chosen parameters.
    """
    sp.run(command.split(" "))


def jellyfish_stats(command: str) -> int:
    """A jellyfish command to get the count which was

    :param command: The jellyfish command with all chosen parameters
    :return: The count of all distinct k-mers.
    """
    result = sp.run(command.split(" "), stdout=sp.PIPE, text=True)
    return int(result.stdout.split("\n")[1].replace(" ", "").split(":")[1])


def count_k_meres(sequence_dict, k=21):
    """Counts all k-meres in the sequences using jellyfish.

    :param sequence_dict: Dictionary with all sequence paths.
    :type sequence_dict: dict[str, str]
    :param k: K-mer length.
    :type k: int
    :return: Species names and number of distinct k-mere.
    """
    k_mer_of_species = list()
    count = 0

    # Iterate through all species.
    for species_name, file_path in sequence_dict.items():
        num_files_to_count = len(sequence_dict) - count
        logger.info("{num} files left to count. Counting {name}", num=num_files_to_count, name=species_name)
        count += 1

        # Set parameters for jellyfish commands.
        k = str(k)
        hash_size = "100M"
        num_threads = "4"
        output_name = str(Path(__file__).parent.parent / "output")

        # Command for jellyfish count.
        count_command = "jellyfish count -m " + k + " -o " + output_name + " -C -s " + hash_size + " -t " + \
                        num_threads + " " + file_path
        # Command for jellyfish stats.
        stats_command = "jellyfish stats " + output_name
        jellyfish_count(count_command)
        k_mer_count = jellyfish_stats(stats_command)

        # Append tuple with species name and distinct k-mer count.
        k_mer_of_species.append((species_name, k_mer_count))

    return k_mer_of_species


def sort_k_mer_counts(k_mer_counts):
    """Sorts the list of k-mers to determine the highest count.

    :param k_mer_counts: List of all species with their k-mer counts.
    :type k_mer_counts: list[tuple[str, int]]
    :return: Sorted list beginning with the highest k-mer count.
    """
    # Define the data type for numpy.
    data_type = [("species", "S50"), ("k_mer_count", int)]
    # Create numpy array with defined data type.
    k_mer_count_sorted = np.array(k_mer_counts, dtype=data_type)
    # Sort array based on k-mer count and than reverse so the first tuple has the highest count.
    k_mer_count_sorted = np.sort(k_mer_count_sorted, order="k_mer_count")[::-1]

    return k_mer_count_sorted


def save_count(k_mer_counts, dir_name):
    """Saves the number of distinct k-meres counted for each species in a txt file.

    :param k_mer_counts: The species with their number of distinct k-meres.
    :type k_mer_counts: list[tuple[str, int]]
    :param dir_name: Name of the parent directory.
    :type dir_name: str
    """
    file_name = Path(__file__).parent.parent / 'genus_metadata' / dir_name / 'distinct_k_meres.txt'
    with open(file_name, "w", encoding='utf-8') as f:
        f.write(f"Species\tdistinct k-meres\n")
        for species in k_mer_counts:
            f.write(f"{species[0].decode('utf-8')}\t{species[1]:,}\n")


def get_highest_k_mer_count(dir_name, k=21):
    """Gets highest k-mer count for all species and k-mer count of genus.

    :param dir_name: Name of the parent directory.
    :type dir_name: str
    :param k: K-mer length.
    :type k: int
    :return: List of the highest k-mer count of all species and the k-mer count of all sequences united.
    """
    # Get highest k-mer count of all species.
    seq_dict = get_seq_paths(dir_name)
    k_mer_counts = count_k_meres(seq_dict, k=k)
    k_mer_sorted = sort_k_mer_counts(k_mer_counts)
    # Uncomment if the k-mer counts should be saved.
    # save_count(k_mer_sorted, dir_name)

    # Count distinct k-mers of genus.
    genus = dir_name.split("_")[0]
    file_name = genus + ".fasta"
    file_path = str(Path(__file__).parent.parent / 'genus_metadata' / dir_name / file_name)
    seq_dict = {genus: file_path}
    k_mer_count = count_k_meres(seq_dict, k=k)

    # Return highest k-mer count of all species and k-mer count of complete genus.
    return [k_mer_sorted[0][1], k_mer_count[0][1]]


def main():
    dir_name = "Listeria_14_12_2022_21-5-13"
    seq_dict = get_seq_paths(dir_name)

    start = perf_counter()

    end = perf_counter()
    print(f"time: {(end-start)/60}")


if __name__ == "__main__":
    main()
