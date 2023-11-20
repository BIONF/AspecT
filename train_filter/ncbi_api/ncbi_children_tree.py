"""This class uses the NCBI Datasets API to get the taxonomy tree of a given Taxon.

The taxonomy tree consists of only the next children to the parent taxon.
The children are only of the next lower rank of the parent taxon.
"""

__author__ = "Berger, Phillip"

import requests

from loguru import logger


class NCBIChildrenTree:
    _taxon: str
    _response: dict
    _parent_taxon_id: str
    _children_taxon_ids = list()
    
    def __init__(self, taxon: str):
        self._taxon = taxon
        self._request_tree()
        
    def _request_tree(self):
        """Make the request for the children tree at the NCBI Datasets API."""
        api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/taxonomy/taxon/{self._taxon}/filtered_subtree"
        raw_response = requests.get(api_url)
        self._response = raw_response.json()["edges"]
        self._parent_taxon_id = str(self._response["1"]["visible_children"][0])
        try:
            tmp_children_ids = self._response[self._parent_taxon_id]["visible_children"]
        except KeyError:
            logger.error("KeyError for key: {key}", key=self._parent_taxon_id)
            logger.error("Response: {response}", response=str(self._response))
            logger.error("Aborting")
            exit()
        for child_id in tmp_children_ids:
            self._children_taxon_ids.append(str(child_id))
    
    def parent_id(self):
        """The NCBI taxon ID of the given Taxon at the initialisation.

        :return: The taxon ID.
        """
        return self._parent_taxon_id

    def children_ids(self) -> list[str]:
        """The NCBI taxon IDs of all children of the given Taxon. The children are all of a lower rank than the parent.

        :return: The taxon IDs as a list.
        """
        return self._children_taxon_ids
        

def main():
    taxon = "286"
    print(NCBIChildrenTree(taxon).children_ids())


if __name__ == "__main__":
    main()
