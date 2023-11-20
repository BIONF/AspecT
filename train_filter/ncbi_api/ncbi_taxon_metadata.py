"""

"""

__author__ = "Berger, Phillip"

import requests

from loguru import logger


class NCBITaxonMetadata:
    _taxon: str
    _response: dict
    _all_metadata: dict

    def __init__(self, taxon: list[str]):
        self._taxon = ",".join(taxon)
        self._all_metadata = dict()
        self._request_metadata()
        self._collect_all_metadata()

    def _request_metadata(self):
        api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/taxonomy/taxon/{str(self._taxon)}"
        raw_response = requests.get(api_url)
        self._response = raw_response.json()["taxonomy_nodes"]

    def _collect_all_metadata(self):
        for child_metadata in self._response:
            taxonomy = child_metadata["taxonomy"]
            rank = taxonomy["rank"]
            name = taxonomy["organism_name"]
            tax_id = str(taxonomy["tax_id"])
            if "Candidatus" not in name:
                if " sp. " not in name:
                    metadata = {
                        "sci_name": name,
                        "tax_id": tax_id,
                        "rank": rank
                    }
                    self._all_metadata[tax_id] = metadata
                else:
                    logger.debug("{name} was not used for training", name=name)
            else:
                logger.debug("{name} was not used for training", name=name)

    def get_response(self):
        return self._response

    def get_metadata(self):
        return self._all_metadata


def main():
    pass


if __name__ == "__main__":
    main()

