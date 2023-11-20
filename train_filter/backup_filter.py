import os
import pickle
from pathlib import Path
from shutil import copy2
from time import asctime, localtime

from loguru import logger

def get_current_time():
    """Returns the current time in the form hh:mm:ss."""
    return asctime(localtime()).split()[3]


def get_dir_names_path():
    """Get the path to the directory_names.txt file which contains the dir_names for the directories with metadata
    for the current bloomfilters.

    :return: Path as string to the directory_names.txt file.
    """
    file_name = "directory_names.txt"

    # Test if file already exists and creates it if not.
    path = Path(__file__).parent.parent / "genus_metadata"
    file_path = path / file_name
    files = os.listdir(path)
    if file_name not in files:
        new_dictionary = {}
        with open(file_path, "wb") as f:
            pickle.dump(new_dictionary, f)

    return str(file_path)


def load_dir_name(genus):
    """Load the dir_names with the metadata for the current bloomfilters. Returns an empty string if the genus does not
    already have bloomfilters.

    :param genus: Name of the genus.
    :type genus: str
    :return: The dir_name where the metadata is for the old bloomfilters of the genus.
    """
    path = get_dir_names_path()

    # Load existing dir_names and add new dir_name.
    dir_names = dict()
    with open(path, "rb") as f:
        dir_names = pickle.load(f)

    if genus in dir_names.keys():
        return dir_names[genus]
    else:
        return ""


def save_dir_name(dir_name):
    """Saves the directory name of the parent directory for the current genus where the assemblies are saved.

    :param dir_name: Name of the directory.
    :type dir_name: str
    """
    genus = dir_name.split("_")[0]
    path = get_dir_names_path()

    # Load existing dir_names and add new dir_name.
    dir_names = dict()
    with open(path, "rb") as f:
        dir_names = pickle.load(f)
    dir_names[genus] = dir_name

    # Save dir_names.
    with open(path, "wb") as f:
        pickle.dump(dir_names, f)


def backup_old_bf(old_dir_name):
    """Backup old bloomfilters so new ones can be created.

    :param old_dir_name: Directory name to save old bloomfilters.
    :type old_dir_name: str
    """
    genus = old_dir_name.split("_")[0]
    # Create backup directory.
    backup_path = Path(__file__).parent.parent

    backup_path = backup_path / "old_filter" / old_dir_name
    backup_path_species = backup_path / genus
    backup_path_sizes = backup_path / "array_sizes"
    backup_path_meta = backup_path / "Metagenomes"
    backup_path_names = backup_path / "species_names"
    backup_path_transl_dict = backup_path / "translation_dicts"

    try:
        os.mkdir(backup_path)
        os.mkdir(backup_path_species)
        os.mkdir(backup_path_sizes)
        os.mkdir(backup_path_meta)
        os.mkdir(backup_path_names)
        os.mkdir(backup_path_transl_dict)
    except FileExistsError:
        logger.info("The backup for the filters of {dir} already exists", dir=old_dir_name)
        logger.info("Stopping backup")
        return None

    # Moving old bloomfilters and extra files.
    # Setting path.
    from_path = Path(__file__).parent.parent / "filter"

    # Copying files.
    # Species filters.
    from_path_species = from_path / genus
    files = os.listdir(from_path_species)
    # Copy all files to back up directory.
    for file in files:
        file_path = from_path_species / str(file)
        copy2(file_path, backup_path_species)

    # Array sizes.
    from_path_sizes = from_path / "array_sizes" / (str(genus) + ".txt")
    copy2(from_path_sizes, backup_path_sizes)

    # Metagenome filter.
    from_path_meta = from_path / "Metagenomes" / (str(genus) + ".txt")
    copy2(from_path_meta, backup_path_meta)

    # Species names.
    from_path_species_names = from_path / "species_names" / ("Filter" + str(genus) + ".txt")
    from_path_meta_names = from_path / "species_names" / ("Filter" + str(genus) + "Complete.txt")
    copy2(from_path_species_names, backup_path_names)
    copy2(from_path_meta_names, backup_path_names)

    # Support vector machine.
    from_path_svm = Path(__file__).parent.parent / "Training_data" / (genus + "_Training_data_spec.csv")
    copy2(from_path_svm, backup_path)

    # Translation dict
    from_path_transl_dict = from_path / "translation_dicts" / (genus + ".pickle")
    # Check if file exits.
    if os.path.exists(from_path_transl_dict):
        copy2(from_path_transl_dict, backup_path_transl_dict)


def main():
    genus = "Salmonella"
    dir_name = "Salmonella_14_1_2023_0-7-23"
    backup_old_bf(dir_name)


if __name__ == "__main__":
    main()
