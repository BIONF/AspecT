"""

"""

__author__ = "Berger, Phillip"

import os

import requests

from pathlib import Path


def download_assemblies(accessions, dir_name, target_folder, zip_file_name):
    """Download and save assemblies from the NCBI database.

    :param accessions: All collected accessions from the NCBI RefSeq-database.
    :type accessions: list
    :param dir_name: Name of the directory where the assemblies will be saved.
    :type dir_name: str
    :param target_folder: Name for the folder in which the downloaded files will be stored.
    :type target_folder: str
    :param zip_file_name: Name of the zip file. E.g. Klebsiella aerogenes.zip.
    :type zip_file_name: str
    """
    path = Path(__file__).parent.parent.parent / 'genus_metadata' / dir_name / target_folder / zip_file_name
    api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/{','.join(accessions)}/download"
    parameters = {
        "include_annotation_type": "GENOME_FASTA",
        "filename": zip_file_name
    }
    os.makedirs(os.path.dirname(path), exist_ok=True)
    genome_download = requests.get(api_url, params=parameters)
    with open(path, "wb") as f:
        f.write(genome_download.content)


def main():
    pass


if __name__ == "__main__":
    main()
