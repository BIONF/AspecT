"""

"""

__author__ = "Berger, Phillip"

import os
from pathlib import Path


def check_folder(folder_path: Path) -> bool:
    """Checks if a folder exists.

    :param folder_path: Path to the folder.
    :return: The result of the test.
    """
    return os.path.isdir(folder_path)


def create_folder(folder_path: Path):
    """Creates a folder with the given path.

    :param folder_path: The path to the folder
    """
    os.mkdir(folder_path)


def check_folder_structure():
    """Checks the folder structure and creates new folders if needed."""
    # Create list of all folder paths.
    path = Path(__file__).parent.parent
    filter_path = path / "filter"
    meta_path = path / "genus_metadata"
    old_filter_path = path / "old_filter"
    filter_folders = ["array_sizes", "Metagenomes", "species_names", "translation_dicts"]
    folder_paths = [filter_path, meta_path, old_filter_path]
    for folder_name in filter_folders:
        filter_folder_path = filter_path / folder_name
        folder_paths.append(filter_folder_path)

    # Check if folders exist. If not create them.
    for folder_path in folder_paths:
        folder_exists = check_folder(folder_path)
        if not folder_exists:
            create_folder(folder_path)


def main():
    pass


if __name__ == "__main__":
    main()
