"""

"""

__author__ = "Berger, Phillip"

import os
import shutil
import zipfile
from pathlib import Path
from time import asctime, localtime

from loguru import logger

class ExtractConcatenate:
    _path = Path(__file__).parent.parent / "genus_metadata"
    _fasta_endings = ["fasta", "fna", "fa", "ffn", "frn"]
    _dir_name: str
    _delete: bool

    def __init__(self, dir_name: str, delete: bool):
        self._path = self._path / dir_name
        self._dir_name = dir_name
        self._delete = delete

    def bf(self):
        zip_path = self._path / "zip_files"
        unzipped_path = self._path / "zip_files_extracted"
        concatenate_path = self._path / "concatenate"
        logger.info("Extracting assemblies")
        self._extract_zip(zip_path, unzipped_path)
        logger.info("Concatenating assemblies")
        self._concatenate_bf(unzipped_path, concatenate_path)
        if self._delete:
            logger.info("Deleting copies")
            self._delete_dir(unzipped_path)

    def svm(self, species_accessions: dict):
        zip_path = self._path / "training_data_zipped"
        unzipped_path = self._path / "training_data_unzipped"
        assemblies_path = self._path / "training_data"
        logger.info("Extracting assemblies")
        self._extract_zip(zip_path, unzipped_path)
        logger.info("Copying assemblies into one folder")
        self._copy_assemblies(unzipped_path, assemblies_path)
        logger.info("Changing headers")
        self._change_header(assemblies_path, species_accessions)
        if self._delete:
            logger.info("Deleting copies")
            self._delete_zip_files(zip_path)
            self._delete_dir(zip_path)
            self._delete_dir(unzipped_path)

    @staticmethod
    def _extract_zip(zip_path, unzipped_path):
        # Make new directory.
        os.mkdir(unzipped_path)

        file_names = os.listdir(zip_path)
        for file in file_names:
            file_path = zip_path / file
            if zipfile.is_zipfile(file_path):
                with zipfile.ZipFile(file_path) as item:
                    directory = unzipped_path / file.replace(".zip", "")
                    item.extractall(directory)

    def _copy_assemblies(self, unzipped_path, assemblies_path):
        os.mkdir(assemblies_path)
        for folder in os.listdir(unzipped_path):
            for (root, dirs, files) in os.walk(unzipped_path / str(folder)):
                for file in files:
                    file_ending = file.split(".")[-1]
                    if file_ending in self._fasta_endings:
                        file_path = Path(root) / file
                        shutil.copy(file_path, (assemblies_path / file))

    @staticmethod
    def _change_header(assemblies_path, species_accessions: dict):
        files = os.listdir(assemblies_path)
        # Iterate through all species.
        for name, accessions in species_accessions.items():
            # Iterate through all accessions of the current species.
            for accession in accessions:
                # Iterate through all file names.
                for file in files:
                    if accession in file:
                        file_path = assemblies_path / str(file)
                        # Change the header.
                        with open(file_path, "r") as f:
                            sequence = ""
                            for line in f.readlines():
                                if line[0] != ">":
                                    sequence += line
                            new_header = f">{name}\n"
                        with open(file_path, "w") as f:
                            f.write(new_header)
                            f.write(sequence)

    def _concatenate_bf(self, unzipped_path, concatenate_path):
        # Make new directory
        os.mkdir(concatenate_path)

        genus = self._dir_name.split("_")[0]
        meta_path = self._path / (genus + ".fasta")
        # Open the meta-genome file.
        with open(meta_path, "w") as meta_file:
            # Write the header.
            meta_header = f">{genus} metagenome\n"
            meta_file.write(meta_header)

            # Open the fasta files for each species.
            for folder in os.listdir(unzipped_path):
                species_files = list()
                # Walk through dirs to get all fasta files.
                for (root, dirs, files) in os.walk(unzipped_path / folder):
                    for file in files:
                        file_ending = file.split(".")[-1]
                        if file_ending in self._fasta_endings:
                            species_files.append(Path(root) / file)

                # Gather all sequences and headers.
                sequences = list()
                headers = list()
                for file in species_files:
                    with open(file, "r") as fasta_file:
                        lines = fasta_file.readlines()
                        header, sequence = self._get_header_and_sequence(lines)
                        headers.append(header)
                        sequences.append(sequence)

                # Concatenate sequences
                species_sequence = "".join(sequences)
                species_header = ">" + " § ".join(headers) + "\n"

                # Save concatenated sequences and headers
                species_path = concatenate_path / (folder + ".fasta")
                with open(species_path, "w") as species_file:
                    species_file.write(species_header)
                    species_file.write(species_sequence)

                # Write the sequence to the meta file.
                meta_file.write(species_sequence)

    @staticmethod
    def _get_header_and_sequence(lines: list):
        sequences = list()
        header = ""
        for line in lines:
            if line[0] == ">":
                header = line[1:].replace("\n", "")
            else:
                line = line.replace("\n", "")
                sequences.append(line)
        sequence = "".join(sequences)
        return header, sequence

    @staticmethod
    def _delete_dir(dir_path):
        shutil.rmtree(dir_path, ignore_errors=False, onerror=None)

    @staticmethod
    def _delete_zip_files(dir_path):
        files = os.listdir(dir_path)
        for file in files:
            if zipfile.is_zipfile(file):
                file_path = dir_path / str(file)
                os.remove(file_path)


def get_current_time():
    """Returns the current time in the form hh:mm:ss."""
    return asctime(localtime()).split()[3]


def main():
    pass


if __name__ == "__main__":
    main()
