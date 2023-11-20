"""

"""

__author__ = "Berger, Phillip"

from time import sleep, asctime, localtime

import requests

from loguru import logger


class NCBIAssemblyMetadata:
    _all_metadata: dict
    _count: int
    _ani_gcf: list
    _parameters: dict
    _accessions: list[str]
    _contig_n50: int
    _all_metadata_complete: dict

    def __init__(self, all_metadata: dict, ani_gcf: list, count=8, contig_n50=10000):
        self._all_metadata = all_metadata
        self._count = count
        self._ani_gcf = ani_gcf
        self._contig_n50 = contig_n50

        self._set_parameters()

        tmp_metadata = dict()
        for tax_id, curr_metadata in self._all_metadata.items():
            sleep(1)
            species_name = curr_metadata['sci_name']
            logger.info("Collecting metadata of {name}", name=species_name)
            accessions = self._make_request(taxon=tax_id)
            if len(accessions) != 0:
                curr_metadata["accessions"] = accessions
                tmp_metadata[tax_id] = curr_metadata

        self._all_metadata_complete = tmp_metadata

    def _set_parameters(self):
        params = {
            "filters.reference_only": "false",
            "filters.assembly_source": "refseq",
            "filters.exclude_atypical": "true",
            "page_size": self._count,
            "page_token": ""
        }
        params_ref = params.copy()
        params_ref["filters.reference_only"] = "true"

        params_comp_genome = params.copy()
        params_comp_genome["filters.assembly_level"] = "complete_genome"

        params_chrom = params.copy()
        params_chrom["filters.assembly_level"] = "chromosome"

        params_scaffold = params.copy()
        params_scaffold["filters.assembly_level"] = "scaffold"

        params_contig = params.copy()
        params_contig["filters.assembly_level"] = "contig"

        self._parameters = {
            "params_ref": params_ref,
            "params_comp_genome": params_comp_genome,
            "params_chrom": params_chrom,
            "params_scaffold": params_scaffold,
            "params_contig": params_contig
        }

    def _make_request(self, taxon: str):
        api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/genome/taxon/{taxon}"
        accessions = list()
        count = 0
        for request_type, parameters in self._parameters.items():
            raw_response = requests.get(api_url, params=parameters)
            response = raw_response.json()
            if response:
                try:
                    assemblies = response["assemblies"]
                    for assembly in assemblies:
                        curr_assembly = assembly["assembly"]
                        curr_accession = curr_assembly["assembly_accession"]
                        curr_contig_n50 = curr_assembly["contig_n50"]
                        if count < self._count:
                            if curr_accession in self._ani_gcf and curr_contig_n50 > self._contig_n50:
                                accessions.append(curr_accession)
                                count += 1
                        else:
                            break
                except KeyError:
                    pass
                    logger.debug("While requesting: {type} an error response was given", type=request_type)
                    logger.debug(str(response))

            if count >= self._count:
                break
        return accessions

    def get_all_metadata(self):
        return self._all_metadata_complete


def get_current_time():
    """Returns the current time in the form hh:mm:ss."""
    return asctime(localtime()).split()[3]


def main():
    pass


if __name__ == "__main__":
    main()
