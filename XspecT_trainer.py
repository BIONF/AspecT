import argparse
import os
import shutil
from linecache import getline
from pathlib import Path
import pickle
from platform import python_version
import sys
from time import localtime, perf_counter, asctime, sleep

from loguru import logger
from numpy import mean


from train_filter.ncbi_api import ncbi_assembly_metadata, ncbi_taxon_metadata, ncbi_children_tree, download_assemblies

from train_filter import create_svm, html_scrap, extract_and_concatenate, get_paths, \
    interface_XspecT, k_mer_count, backup_filter, check_folders


def check_user_input(user_input: str):
    """The given input of the user will be checked. The input has to be a genus in NCBI.

    :return: The genus name.
    """
    taxon_metadata = ncbi_taxon_metadata.NCBITaxonMetadata([user_input])
    all_metadata = taxon_metadata.get_metadata()
    for metadata in all_metadata.values():
        sci_name = metadata["sci_name"]
        tax_id = metadata["tax_id"]
        rank = metadata["rank"]
        if not sci_name == user_input and not tax_id == user_input:
            print(f"{get_current_time()}| The given genus: {user_input} was found as genus: {sci_name} "
                  f"ID: {tax_id}")
            print(f"{get_current_time()}| Using {sci_name} as genus name.")
        if rank == "GENUS":
            return str(sci_name)
        else:
            print(f"{get_current_time()}| {user_input} is rank {rank} and not genus.")
            exit()


def copy_custom_data(bf_path: str, svm_path: str, dir_name: str):
    """

    :param bf_path:
    :param svm_path:
    :param dir_name:
    :return:
    """
    path = Path(__file__).parent.absolute() / "genus_metadata" / dir_name
    new_bf_path = path / "concatenate"
    new_svm_path = path / "training_data"

    # Make the new directories.
    os.mkdir(path)
    os.mkdir(new_bf_path)
    os.mkdir(new_svm_path)

    # Move bloomfilter files.
    bf_files = os.listdir(bf_path)
    for file in bf_files:
        file_path = Path(bf_path) / file
        new_file_path = new_bf_path / file
        shutil.copy2(file_path, new_file_path)

    # Move svm files.
    svm_files = os.listdir(svm_path)
    for file in svm_files:
        file_path = Path(svm_path) / file
        new_file_path = new_svm_path / file
        shutil.copy2(file_path, new_file_path)


def concatenate_meta(dir_name: str):
    """ Concatenates all concatenated fasta files that are used to train bloomfilters to one fasta file.

    :param dir_name:
    :return:
    """
    path = Path(__file__).parent.absolute() / "genus_metadata" / dir_name
    files_path = path / "concatenate"
    fasta_endings = ["fasta", "fna", "fa", "ffn", "frn"]
    genus = dir_name.split("_")[0]
    meta_path = path / (genus + ".fasta")
    files = os.listdir(files_path)

    with open(meta_path, "w") as meta_file:
        # Write the header.
        meta_header = f">{genus} metagenome\n"
        meta_file.write(meta_header)

        # Open each concatenated species file and write the sequence in the meta file.
        for file in files:
            file_ending = str(file).split(".")[-1]
            if file_ending in fasta_endings:
                with open((files_path / str(file)), "r") as species_file:
                    for line in species_file:
                        if line[0] != ">":
                            meta_file.write(line.replace("\n", ""))


def count_avg_seq_len(dir_name):
    """Counts the sequence length for each species concatenated fasta file and computes the average.

    :param dir_name: Directory name for current genus.
    :type dir_name: str
    :return: The average sequence length.
    """
    path = get_paths.get_concatenate_file_path(dir_name)

    # Create list with all sequence lengths.
    files = os.listdir(path)
    counts = list()
    for file in files:
        file_path = path / str(file)
        sequence = getline(str(file_path), 2)
        counts.append(len(sequence))

    # Return avg. sequence length.
    return int(round(float(mean(counts)), 0))


def check_meta_file_size(dir_name) -> bool:
    """ Checks the metagenome fasta file if every concatenated species file was used for
    its creation by comparing the file sizes.

    :param dir_name: Directory name for current genus.
    :type dir_name: str
    :return: True or False depending on the answer.
    """
    path = Path(__file__).parent.absolute()
    species_path = path / "genus_metadata" / dir_name / "concatenate"
    genus = dir_name.split("_")[0]
    meta_path = path / "genus_metadata" / dir_name / (str(genus) + ".fasta")

    species_files = os.listdir(species_path)
    all_files_size = 0
    for file in species_files:
        file_size = os.path.getsize(species_path / str(file))
        all_files_size += file_size

    meta_file_size = os.path.getsize(meta_path)

    all_files_size = round(all_files_size / (1024 ** 2))
    meta_file_size = round(meta_file_size / (1024 ** 2))

    # Compare both sizes.
    same = False
    if all_files_size == meta_file_size:
        same = True
    return same


def check_meta_file_content(dir_name: str):
    """Checks if every sequence used to concatenate the meta file is fully inside the meta file.

    :param dir_name: Directory name for current genus.
    :return: True if every sequence is inside the meta file.
    """
    path = Path(__file__).parent / "genus_metadata" / dir_name
    concatenate_path = path / "concatenate"
    genus = dir_name.split("_")[0]
    mg_file_name = f"{genus}.fasta"
    mg_file_path = path / mg_file_name
    mg_str = ""
    with open(mg_file_path, "r") as mg_file:
        for line in mg_file:
            if line[0] != ">":
                mg_str = line
    files = os.listdir(concatenate_path)

    for file in files:
        file_path = concatenate_path / file
        with open(file_path, "r") as con_file:
            for line in con_file:
                if line[0] == ">":
                    continue

                if line not in mg_str:
                    logger.error(f"{file} not in metagenome")
                else:
                    logger.info(f"{file} in metagenome")


def init_argparse() -> argparse.ArgumentParser:
    """Initiate the command line parser for XspecT-trainer.py"""
    parser = argparse.ArgumentParser(
        prog="XspecT-trainer",
        description="Automatically trains bloomfilter, of a given genus, so they can later be used by XspecT to "
                    "assign species to assemblies."
    )
    parser.add_argument("genus", type=str,
                        help="The name of the genus for which the filters will be trained. Can also be a NCBI Taxon ID")
    parser.add_argument("mode", metavar="mode", choices=["1", "2", "3"], type=str,
                        help="Declares which mode should be used. 1: Train bloomfilters with assemblies from the ncbi"
                             "RefSeq database. 2: Train bloomfilters with custom assembles. The paths to the assemblies"
                             "need to be given. 3: Check if metagenome file was correctly created.")
    parser.add_argument("-bf", "--bf_path", type=str,
                        help="The path to the assemblies that will be used to train the bloomfilters if mode 2 is used."
                        )
    parser.add_argument("-svm", "--svm_path", type=str,
                        help="The path to the assemblies that will be used to train the support-vector-machine if "
                             "mode 2 is used.")
    parser.add_argument("-c", "--complete", action="store_true",
                        help="Declares if all of every 500th k-mere should be used to train the bloomfilters.")
    parser.add_argument("-d", "--dir_name", type=str,
                        help="Write the directory name from genus_metadata to check metagenome file.")

    return parser


def set_logger(dir_name: str):
    """Sets the logger parameters.

    :param dir_name: Name of the folder where the log should be saved.
    """
    genus = dir_name.split("_")[0]

    # Starting logger.
    logger.remove()
    logger.add(sys.stderr, format="{time:HH:mm:ss} | {level} | {message}", level="INFO")
    log_path = (Path(__file__).parent / "genus_metadata" / dir_name / (genus + ".log"))
    logger.add(log_path, format="{time:HH:mm:ss} | {level} | {message}", level="DEBUG")


def create_translation_dict(dir_name: str) -> dict[str, str]:
    """Create a translation dictionary to translate the taxon ID to its scientific name from the file names.

    :param dir_name: Directory name for current genus.
    :return: The created translation dictionary.
    """
    path = Path(__file__).parent / "genus_metadata" / dir_name / "concatenate"
    files = os.listdir(path)
    translation_dict = dict()
    for file in files:
        file_split = file.split(".")[0].split("_")
        tax_id = file_split[0]
        name = file_split[1]
        translation_dict[tax_id] = name

    return translation_dict


def save_translation_dict(dir_name: str, translation_dict: dict[str, str]):
    """Saves the translation dict in filter/translation_dicts as pickle file.

    :param dir_name: Directory name for current genus.
    :param translation_dict: A dictionary with taxon ID as key and its corresponding scientific name as value.
    """
    genus = dir_name.split("_")[0]
    folder_path = Path(__file__).parent / "filter" / "translation_dicts"
    # Check if folder exists
    if os.path.exists(folder_path):
        # Check if it is a folder
        if not os.path.isdir(folder_path):
            logger.error("Path: {path} is not a folder", path=folder_path)
            logger.error("Aborting")
            exit()
    else:
        # Create folder
        os.mkdir(folder_path)

    file_name = f"{genus}.pickle"
    file_path = folder_path / file_name
    with open(file_path, "wb") as f:
        pickle.dump(translation_dict, f)


def change_bf_assembly_file_names(dir_name: str):
    """ Change all concatenated assembly names to only the taxon ID.

    :param dir_name: Directory name for current genus.
    """
    path = Path(__file__).parent / "genus_metadata" / dir_name / "concatenate"
    files = os.listdir(path)
    for file in files:
        file_split = file.split(".")[0].split("_")
        tax_id = file_split[0]
        new_file_name = f"{tax_id}.fasta"
        os.rename((path / file), (path / new_file_name))


def get_current_time():
    """Returns the current time in the form hh:mm:ss."""
    return asctime(localtime()).split()[3]


def main():
    # command line should look like this: python XspecT_trainer.py genus mode path_to_bf_files path_to_svm_files
    parser = init_argparse()
    args = parser.parse_args()
    genus = args.genus
    mode = args.mode
    complete = args.complete
    if complete:
        spacing = 1
    else:
        spacing = 500

    # Check folder structure
    check_folders.check_folder_structure()

    # Check user input.
    genus = check_user_input(user_input=genus)

    # The directory name is defined in the following format: 'genus'_DD_MM_YYYY_hh-mm-ss
    curr_time = localtime()
    dir_name = f"{genus}_{curr_time[2]}_{curr_time[1]}_{curr_time[0]}_{curr_time[3]}-{curr_time[4]}-{curr_time[5]}"

    # Set the logger.
    set_logger(dir_name)

    # Time for the whole program.
    start_all = perf_counter()
    if mode == "1":
        # Search for every defined species of the genus.
        start_tax = perf_counter()
        logger.info("Getting all species of the genus")
        children_ids = ncbi_children_tree.NCBIChildrenTree(genus).children_ids()
        species_dict = ncbi_taxon_metadata.NCBITaxonMetadata(children_ids).get_metadata()
        end_tax = perf_counter()

        # Get all gcf accessions that have Taxonomy check result OK.
        logger.info("Checking ANI data for updates")
        ani = html_scrap.TaxonomyCheck()
        ani_gcf = ani.ani_gcf()

        # Look for up to 8 assembly accessions per species.
        start_meta = perf_counter()
        logger.info("Getting assembly metadata")
        all_metadata = ncbi_assembly_metadata.NCBIAssemblyMetadata(
            all_metadata=species_dict,
            ani_gcf=ani_gcf,
            count=8,
            contig_n50=10000
        )
        all_metadata = all_metadata.get_all_metadata()
        logger.info("Finished metadata collecting\n")
        end_meta = perf_counter()

        # Download the chosen assemblies.
        # One file for each species with it's downloaded assemblies in zip format.
        start_download = perf_counter()
        """for key, value in all_metadata.items():
            print(f"{key}: {value}")"""

        # Iterate through all species.
        logger.info("Downloading assemblies for bloomfilter training")
        for metadata in all_metadata.values():
            # Only try to download when the species has accessions.
            if len(metadata["accessions"]) >= 1:
                sleep(5)
                species_name = metadata["sci_name"]
                tax_id = metadata['tax_id']
                logger.info("Downloading {id}_{name}", id=tax_id, name=species_name)
                file_name = f"{tax_id}_{species_name}.zip"

                # Selecting the first 4 assemblies for training the filters.
                accessions = list()
                for accession in metadata["accessions"]:
                    accessions.append(accession)
                    if len(accessions) == 4:
                        break

                download_assemblies.download_assemblies(
                    accessions=accessions,
                    dir_name=dir_name,
                    target_folder="zip_files",
                    zip_file_name=file_name
                )
        logger.info("Downloads finished\n")
        end_download = perf_counter()

        # Concatenate all assemblies of each species.
        start_concatenate = perf_counter()
        extract_bf = extract_and_concatenate.ExtractConcatenate(
            dir_name=dir_name,
            delete=True
        )
        extract_bf.bf()
        logger.info("Finished extracting and concatenating\n")
        end_concatenate = perf_counter()

        # Compute average sequence length.
        avg_len = count_avg_seq_len(dir_name)

        # Download assemblies for svm creation.
        start_svm_dl = perf_counter()
        logger.info("Downloading assemblies for support-vector-machine training")
        accessions = dict()
        for metadata in all_metadata.values():
            # Only add taxon with accessions.
            if len(metadata["accessions"]) >= 1:
                accessions[metadata["tax_id"]] = metadata["accessions"]

        # Downloading assemblies.
        create_svm.get_svm_assemblies(
            all_accessions=accessions,
            dir_name=dir_name
        )
        logger.info("Finished downloading\n")

        # Extracting assemblies.
        extract_svm = extract_and_concatenate.ExtractConcatenate(
            dir_name=dir_name,
            delete=True
        )
        extract_svm.svm(species_accessions=accessions)

        end_svm_dl = perf_counter()

    elif mode == "2":
        # Mode 2 needs to folders one with concatenated fasta files.
        # The files should have .fasta as a file ending and its name should be the species ID and its name without
        # the genus name. e.g. 28901_enterica.fasta for Salmonella enterica. The ID can be any ID. The standard is ncbi
        # taxon IDs. The ID should only contain numbers from 0-9.
        # The second folder should have assembly fasta files for every species. These should have a code in the file
        # name to understand where the data came from. Its header should have > at the start and after the species
        # ID. E.g. >28901\n
        # Path to concatenated fasta files for bloomfilter training.
        bf_path = args.bf_path
        # Path to fasta files for svm training.
        svm_path = args.svm_path

        # Check if paths were given.
        if bf_path:
            if not os.path.exists(bf_path):
                logger.error("The given path to the bloomfilter assemblies doesn't exist")
                logger.error("Aborting")
                exit()
        else:
            logger.error("There was no path to the bloomfilter assemblies given")
            logger.error("Aborting")
            exit()
        if svm_path:
            if not os.path.exists(svm_path):
                logger.error("The given path to the support-vector-machine assemblies doesn't exist")
                logger.error("Aborting")
                exit()
        else:
            logger.error("There was no path to the support-vector-machine assemblies given")
            logger.error("Aborting")
            exit()

        # Move the given files to genus_metadata.
        logger.info("Copying data given into genus_metadata")
        copy_custom_data(bf_path=bf_path,
                         svm_path=svm_path,
                         dir_name=dir_name)

        # Create Metagenome fasta file of all concatenated fasta files.
        logger.info("Creating meta fasta file")
        concatenate_meta(dir_name)

    elif mode == "3":
        logger.info("Checking metagenome file")
        mg_check_dir_name = args.dir_name
        if not mg_check_dir_name:
            logger.error("There was no directory name given")
            logger.error("Aborting")
            exit()
        check_meta_file_content(mg_check_dir_name)
        logger.info("Finished")
        logger.opt(record=True).info("Elapsed time: {record[elapsed]}")
        exit()

    # Check file sizes.
    result = False
    logger.info("Checking if metagenome file was correctly created")
    result = check_meta_file_size(dir_name)
    count = 0
    while not result:
        logger.error("Metagenome file was not correctly created")
        logger.info("Trying to remake metagenome fasta file")
        concatenate_meta(dir_name)
        logger.info("Rechecking metagenome file")
        result = check_meta_file_size(dir_name)
        count += 1
        if count == 3:
            logger.error("Can't create metagenome file")
            logger.error("Aborting")
            exit()

    # Make dictionary for translating taxon ID to scientific name.
    translation_dict = create_translation_dict(dir_name)
    change_bf_assembly_file_names(dir_name)

    # Count all distinct k-mers and return the highest count.
    start_count = perf_counter()
    logger.info("Counting all distinct k-meres")
    highest_counts = k_mer_count.get_highest_k_mer_count(dir_name)
    end_count = perf_counter()

    # Backup of old bloomfilters if genus already has them.
    old_dir_name = backup_filter.load_dir_name(genus)
    # Only backup if old bloomfilters exist.
    logger.info("Doing old filter backup")
    if old_dir_name:
        backup_filter.backup_old_bf(old_dir_name)
    else:
        logger.info("No backup was needed")

    # Train new Bloomfilters with concatenated files of each species.
    start_bf = perf_counter()
    # Compute the array size with the highest count of distinct k-mers.
    array_size_species = int(round(interface_XspecT.compute_array_size(highest_counts[0]) + 1000000, -6))
    array_size_complete = int(round(interface_XspecT.compute_array_size(highest_counts[1]) + 1000000, -6))

    # Save array sizes for XspecT.
    logger.info("Saving bloomfilter sizes\n")
    interface_XspecT.save_array_sizes(genus, [str(array_size_species), str(array_size_complete)])

    # Train Bloomfilters of species.
    logger.info("Training bloomfilters")
    files_path, result_path = interface_XspecT.make_paths(dir_name, genus)
    interface_XspecT.new_train_core(files_path, result_path, array_size_species)
    interface_XspecT.new_write_file_dyn(result_path, genus, meta_mode=False)

    # Train Bloomfilter for complete genus.
    logger.info("Training metagenome bloomfilter")
    files_path = get_paths.get_current_dir_file_path(dir_name)

    result_path = get_paths.get_metagenome_filter_path()

    interface_XspecT.new_train_core(str(files_path), str(result_path), array_size_complete)
    interface_XspecT.new_write_file_dyn(str(result_path), genus, meta_mode=True)

    end_bf = perf_counter()

    # Create support vector machine.
    start_svm = perf_counter()
    logger.info("Training support-vector-machine")
    # Create svm.
    create_svm.new_helper(
        spacing,
        genus=genus,
        dir_name=dir_name,
        array_size=array_size_species,
        k=21
    )
    end_svm = perf_counter()
    end_all = perf_counter()

    logger.info("Program runtime: {time} m", time=(round((end_all-start_all)/60, 2)))

    if mode == "1":
        # Print and save collected statistics.
        logger.info("Saving collected runtime statistics")
        time_print = f"Python version: {python_version()} \n" \
                     f"Average sequence length: {avg_len:,} \n" \
                     f"All time: {(end_all-start_all)/60:.2f} m\n" \
                     f"Tax time: {(end_tax-start_tax):.2f} s\n" \
                     f"Meta time: {(end_meta-start_meta)/60:.2f} m\n" \
                     f"Download time: {(end_download-start_download)/60:.2f} m\n" \
                     f"Concatenate time: {(end_concatenate-start_concatenate):.2f} s\n" \
                     f"Count time: {(end_count-start_count)/60:.2f} m\n" \
                     f"Training time: {(end_bf-start_bf)/60:.2f} m\n" \
                     f"Support vector machine time: {((end_svm+end_svm_dl)-(start_svm+start_svm_dl))/60:.2f} m\n"

        # Save time measurements.
        interface_XspecT.save_time_stats(time_print, dir_name)

    # Save translation dict
    save_translation_dict(dir_name, translation_dict)

    # Save directory name for later filter backup.
    backup_filter.save_dir_name(dir_name)

    logger.info("XspecT-trainer is finished.")


if __name__ == "__main__":
    main()
