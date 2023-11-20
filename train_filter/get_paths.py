from pathlib import Path


def get_zip_files_path(dir_name):
    """Returns str of file path to the zip files.

    :param dir_name: Name of the current genus_metadata directory.
    :type dir_name: str
    :return: File path to zip files.
    """
    return Path(__file__).parent.parent / 'genus_metadata' / dir_name / 'zip_files_extracted'


def get_concatenate_file_path(dir_name):
    """Returns str to file path of the concatenate directory.

    :param dir_name: Name of the current genus_metadata directory.
    :type dir_name: str
    :return: File path to the concatenated species assemblies.
    """
    return Path(__file__).parent.parent / 'genus_metadata' / dir_name / 'concatenate'


def get_current_dir_file_path(dir_name):
    """Returns str of file path to the directory with the currently needed metagenome assembly.

    :param dir_name: Name of the current genus_metadata directory.
    :type dir_name: str
    :return: File path to the metagenome assembly.
    """
    return Path(__file__).parent.parent / 'genus_metadata' / dir_name


def get_metagenome_filter_path():
    """Returns the file path to the metagenome filters."""
    return Path(__file__).parent.parent / 'filter' / 'Metagenomes'


def main():
    pass


if __name__ == "__main__":
    main()
