import os
import pickle
from pathlib import Path
from shutil import rmtree

from loguru import logger
from numpy import log, square

import BF_v2


def compute_array_size(n, p=0.01):
    """Computes the Bit-Array-Size for the bloomfilters.

    :param n: Highest k-mer count of a species.
    :type n: int
    :param p: Rate of mistakes.
    :type p: float
    :return: Bit-Array-Size for the bloomfilters.
    """
    return -((n * log(p)) / (square(log(2))))


def make_paths(dir_name, genus):
    """ Create paths to the concatenated sequences and to where the new bloomfilters will be saved.

    :param dir_name: Name of the parent directory.
    :type dir_name: str
    :param genus: Name of the genus.
    :type genus: str
    :return: The path to the sequence files and the bloomfilter directory.
    """
    # Path to concatenated sequences
    files_path = Path(__file__).parent.parent / 'genus_metadata' / dir_name / 'concatenate'

    # Path for results.
    result_path = Path(__file__).parent.parent / 'filter' / genus
    # Try to create the directory for the bloomfilters.
    try:
        os.mkdir(result_path)
    except FileExistsError:
        # Delete the old directory with bloomfilters if already existed.
        rmtree(result_path, ignore_errors=False, onerror=None)
        os.mkdir(result_path)

    return str(files_path), str(result_path)


def init_bf(array_size, clonetypes=1, hashes=7, k=21):
    """Initiates an bloomfilter object  with given parameters.

    :param array_size: The size for the byte-array.
    :type array_size: int
    :param clonetypes: Number of clonetypes.
    :type clonetypes: int
    :param hashes: Number of hash functions used.
    :type hashes: int
    :param k: Length of k-mers.
    :type k: int
    :return: The initiated bloomfilter object.
    """
    BF = BF_v2.AbaumanniiBloomfilter(array_size)
    BF.set_arraysize(array_size)
    BF.set_clonetypes(clonetypes)
    BF.set_hashes(hashes)
    BF.set_k(k)
    return BF


def new_train_core(files_path, result_path, array_size, k=21):
    """Trains concatenated genomes into Bloomfilter and saves them.

    :param files_path: Path to where the concatenated sequences are stored.
    :type files_path: str
    :param result_path: Path where the generated Bloomfilter will be saved.
    :type result_path: str
    :param array_size: Array-size for the Bloomfilter.
    :type array_size: int
    :param k: Length of substring.
    :type k: int
    """
    files = os.listdir(files_path)
    # Iterate the files backwards to delete all non fasta files from the list.
    for i in range(len(files) - 1, -1, -1):
        if 'fna' in files[i] or 'fasta' in files[i]:
            continue
        else:
            del files[i]

    # Train a bloomfilter for each species.
    for i in range(len(files)):
        BF = init_bf(array_size=array_size, clonetypes=1, hashes=7, k=k)
        path = Path(files_path) / files[i]
        species_name = files[i].split('.')[0]
        file_name = species_name + '.txt'
        logger.info("Training {name}", name=species_name)
        result = Path(result_path) / file_name
        BF.train_sequence(path, 0)
        BF.save_clonetypes(result)
        BF.cleanup()


def new_write_file_dyn(bf_path, genus, meta_mode=False):
    """Write file with pickled list of all names for the bloomfilters.

    :param bf_path: Path to the bloomfilters.
    :type bf_path: str
    :param genus: Name of the genus.
    :type genus: str
    :param meta_mode: Declare to which bloomfilters the path leads.
    :type meta_mode: bool
    """
    files = os.listdir(bf_path)
    # If the Bloomfilter path leads to Bloomfilter for the metagenome mode.
    if meta_mode:
        for i in range(len(files) -1, -1, -1):
            if genus not in files[i]:
                del files[i]
            else:
                files[i] = files[i][:-4]
        file_name = 'Filter' + genus + 'Complete.txt'

    # If the path leads to bloomfilters for the species.
    else:
        for i in range(len(files) -1, -1, -1):
            if 'txt' not in files[i]:
                del files[i]
            else:
                files[i] = files[i][:-4]
        file_name = 'Filter' + genus + '.txt'

    # Make path for the txt file.
    file_path = Path(__file__).parent.parent / 'filter' / 'species_names' / file_name
    with open(file_path, 'wb') as fp:
        pickle.dump(sorted(files), fp)


def save_array_sizes(genus, array_sizes):
    """Saves the array sizes of the bytearray for the bloomfilters in a txt file.

    :param genus: The current genus.
    :type genus: str
    :param array_sizes: List of all computed array sizes for this genus.
    :type array_sizes: list[str]
    """
    file_name = genus + ".txt"
    path = Path(__file__).parent.parent / 'filter' / 'array_sizes' / file_name

    # Save both array sizes as a string in the format: 'size1 size2' as a txt file.
    # The first size is of the species level filters and the second of the meta-mode filter.
    text = " ".join(array_sizes)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def save_name_dict(genus, name_dict: dict):
    """Saves the names and taxon IDs of all species for which filter were trained. XspecT uses this dict to switch
    between the species names and it's ID. The dict is saved as a csv file.

    :param genus: The genus for which filters were trained.
    :param name_dict: The dictionary with all species names and taxon IDs
    """



def save_time_stats(time_stats, dir_name):
    """Saves the collected time measurements as a txt file.

    :param time_stats: The collected time measurements as a formatted string.
    :type time_stats: str
    :param dir_name: Name of the parent directory.
    :type dir_name: str
    """
    time_file = Path(__file__).parent.parent / 'genus_metadata' / dir_name / "time.txt"
    with open(str(time_file), "w+", encoding="utf-8") as f:
        f.write(time_stats)


def load_translation_dict(genus: str) -> dict[str, str]:
    """Loads the translation dict for the given genus. The key is the taxon ID and its value the scientific name.

    :param genus: The name of the genus.
    :return: The translation dict for the genus.
    """
    file_name = f"{genus}.pickle"
    path = Path(__file__).parent.parent / "filter" / "translation_dicts" / file_name
    with open(path, "rb") as f:
        translation_dict = pickle.load(f)

    return translation_dict


def main():
    a = 28858023
    b = compute_array_size(a)
    print(int(round(b+1000000, -6)))
    # genera = get_genera_array_sizes()
    # print(f"Species: ")
    # print(pre_process_all(genera, meta_mode=False))
    # print(f"\nMeta: ")
    # print(pre_process_all(genera, meta_mode=True))


if __name__ == "__main__":
    main()
