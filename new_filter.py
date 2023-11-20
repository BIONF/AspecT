import os
import pickle
import shutil
from pathlib import Path


def change_dirs():
    # Filter directories.
    path = Path(__file__).parent.absolute()
    filter_path = path / "filter"
    os.mkdir((filter_path / "array_sizes"))
    os.mkdir((filter_path / "Metagenomes"))
    os.mkdir((filter_path / "species_names"))

    # Meta directory.
    meta_path = path / "genus_metadata"
    os.mkdir(meta_path)

    # Backup directory.
    backup_path = path / "old_filter"
    os.mkdir(backup_path)


def move_files():
    path = Path(__file__).parent.absolute()
    filter_path = path / "filter"

    # Metagenome filter.
    old_meta_path = filter_path / "acinetobacter" / "Acinetobacter_Masterv2.txt"
    new_meta_path = filter_path / "Metagenomes" / "Acinetobacter.txt"
    shutil.move(old_meta_path, new_meta_path)
    shutil.rmtree((filter_path / "acinetobacter"))

    # species filter
    os.rename((filter_path / "species"), (filter_path / "Acinetobacter"))

    # Acinetobacter array sizes.
    array_sizes = "173000000 3080000000"
    sizes_path = filter_path / "array_sizes" / "Acinetobacter.txt"
    with open(sizes_path, "w") as f:
        f.write(array_sizes)

    # Species names.
    old_species_names_path = filter_path / "FilterSpecies.txt"
    new_species_names_path = filter_path / "species_names" / "FilterAcinetobacter.txt"
    shutil.move(old_species_names_path, new_species_names_path)

    # Meta name.
    old_meta_names_path = filter_path / "FilterAcinetobacter.txt"
    os.remove(old_meta_names_path)
    text = "Acinetobacter"
    new_meta_names_path = filter_path / "species_names" / "FilterAcinetobacterComplete.txt"
    with open(new_meta_names_path, "wb") as f:
        pickle.dump(text, f)

    # svm
    svm_path = path / "Training_data"
    old_svm_name = svm_path / "Training_data_spec.csv"
    new_svm_name = svm_path / "Acinetobacter_Training_data_spec.csv"
    os.rename(old_svm_name, new_svm_name)


def main():
    change_dirs()
    move_files()


if __name__ == "__main__":
    main()
