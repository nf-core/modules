#!/usr/bin/env python3
import requests
import sys


class Arguments:
    """
    This class contains the arguments of the module.
    """
    def __init__(self) -> None:
        self.input = "$gene_list"
        self.prefix = (
            "$task.ext.prefix"
            if "$task.ext.prefix" != "null"
            else "$meta.id"
        )
        self.output = self.prefix + ".tsv"

class Handler:
    """
    This class some functions to load and save files.
    """

    @staticmethod
    def load_list(filename: str) -> list:
        """
        It reads the a list from a file.
        Each row is an element of the list.
        """
        with open(filename, "r") as f:
            mylist = f.readlines()
        mylist = [element.strip() for element in mylist]
        return mylist

    @staticmethod
    def save_dictionary_to_file(results: dict, filename: str) -> None:
        """
        It saves the a dictionary of list to a file.
        Each list should have the same length
        """
        # check that all lists have the same length
        keys = list(results.keys())
        num_records = len(results[keys[0]])
        for key in keys:
            if len(results[key]) != num_records:
                raise ValueError("All lists should have the same length")
            
        # save the dictionary to a file
        with open(filename, 'w') as f:
            header = "\\t".join(keys) + "\\n"
            f.write(header)
            for i in range(num_records):
                record = "\\t".join(str(results[key][i]) for key in keys) + "\\n"
                f.write(record)

    @staticmethod
    def get_versions(modules: list) -> dict:
        """
        This function takes a list of modules and returns a dictionary with the versions of each module.
        This dictionary will also contain the version of python.
        """
        dic = {
            **{"python": '.'.join(map(str, sys.version_info[:3]))},
            **{module.__name__: module.__version__ for module in modules}
        }
        return dic
    
    @staticmethod
    def format_yaml_like(data: dict, indent: int = 0) -> str:
        """Formats a dictionary to a YAML-like string.

        Args:
            data (dict): The dictionary to format.
            indent (int): The current indentation level.

        Returns:
            str: A string formatted as YAML.
        """
        yaml_str = ""
        for key, value in data.items():
            spaces = "  " * indent
            if isinstance(value, dict):
                yaml_str += f"{spaces}{key}:\\n{Handler.format_yaml_like(value, indent + 1)}"
            else:
                yaml_str += f"{spaces}{key}: {value}\\n"
        return yaml_str


class QuickGO:
    def __init__(self, gene_id: str) -> None:
        self.gene_id = gene_id
        self.results = self.request_quickgo()
        self.results = self.parse_quickgo()
    
    def request_quickgo(self) -> list:
        """
        This function takes a gene id and requests the information from quickgo server.
        It returns a list of dictionaries. Each dictionary has keys like "geneProductId", "goId", "evidence", "goAspect", etc.
        If the request is not succesful, it raises an error and stops the program.
        """
        requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=" + self.gene_id
        r = requests.get(requestURL, headers={ "Accept" : "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        return r.json()['results']

    def parse_quickgo(self) -> dict:
        """
        This function parses the requested object from QuickGo into a dictionary with the following keys:
        - geneId
        - goIid
        - goEvidence
        - goAspect
        """
        results = {'geneId': [], 'goId': [], 'goEvidence': [], 'goAspect': []}
        for result in self.results:
            results['geneId'].append(self.gene_id)
            results['goId'].append(result['goId'])
            results['goEvidence'].append(result['goEvidence'])
            results['goAspect'].append(result['goAspect'])
        return results


if __name__ == "__main__":
    args = Arguments()

    # read the gene ids from the file
    gene_ids = Handler.load_list(args.input)

    # get the GO terms for each gene id
    results = {"geneId": [], "goId": [], "goEvidence": [], "goAspect": []}
    for gene_id in gene_ids:
        go = QuickGO(gene_id)
        for key,values in go.results.items():
            results[key].extend(values)

    # save results to output
    Handler.save_dictionary_to_file(results, args.output)

    # write versions to file
    versions_this_module = {}
    versions_this_module["${task.process}"] = Handler.get_versions([requests])
    with open("versions.yml", "w") as f:
        f.write(Handler.format_yaml_like(versions_this_module))
