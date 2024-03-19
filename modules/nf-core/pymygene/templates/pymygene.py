#!/usr/bin/env python3
import argparse
import mygene
import shlex
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
        self.output_tsv = self.prefix + ".tsv"
        self.output_gmt = self.prefix + ".gmt"
        self.parse_ext_args("$task.ext.args")

    def parse_ext_args(self, args_string: str) -> None:
        """
        It parses the extended arguments.
        """
        # skip when there are no extended arguments
        if args_string == "null":
            args_string = ""
        
        # Parse the extended arguments
        args_list = shlex.split(args_string)  # Split the string into a list of arguments
        parser = argparse.ArgumentParser()
        parser.add_argument('--species', default=None, help="Comma separated of common name of the species or taxon ids")
        parser.add_argument('--filter_go_category', default=None, help="Comma separated list of GO categories to keep. Default: all")
        parser.add_argument('--filter_go_evidence', default=None, help="Comma separated list of GO evidence codes to keep. Default: all")
        args = parser.parse_args(args_list)

        # Convert "null" values to default values
        for attr in vars(args):
            value = getattr(args, attr)
            if value == "null":
                setattr(args, attr, parser.get_default(attr))

        # Assign args attributes to self attributes
        for attr in vars(args):
            setattr(self, attr, getattr(args, attr))


class Version:
    """
    This class contains the functions to get the versions of the modules used in the script.
    """
    
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
                yaml_str += f"{spaces}{key}:\\n{Version.format_yaml_like(value, indent + 1)}"
            else:
                yaml_str += f"{spaces}{key}: {value}\\n"
        return yaml_str


class MyGene:
    """
    This class contains the functions to query the MyGene.info API.
    """
    def __init__(self, filename, species:str = None, filter_go_category: str = None, filter_go_evidence: str = None) -> None:
        self.query = self.load_list(filename)
        self.species = species
        self.filter_go_category = filter_go_category
        self.filter_go_evidence = filter_go_evidence
        self.mg = mygene.MyGeneInfo()
        self.idmap = self.query2idmap()
        self.gene_based_info = self.parse_gene_based_info()
        self.go_based_info = self.parse_go_based_info()

    def load_list(self, filename: str) -> list:
        """
        It reads the a list from a file.
        Each row is an element of the list.
        """
        with open(filename, "r") as f:
            mylist = f.readlines()
        mylist = [element.strip() for element in mylist]
        return mylist

    def query2idmap(self) -> dict:
        """
        It returns a dictionary with the query ids as keys and the mygene ids as values.
        """
        return {dic['_id']: dic['query'] for dic in self.mg.querymany(self.query, species=self.species) if '_id' in dic}

    def parse_gene_based_info(self) -> dict:
        """
        It returns the parsed information for all the query ids.
        The returned information is a dictionary {query gene: {}} of dictionaries with the following keys:
            - query
            - mygene_id
            - go_id
            - go_term
            - go_evidence
            - go_category
            - symbol
            - name
            - taxid
        Note that this is a query gene centric dictionary.
        """
        info = {}
        for dic in self.mg.getgenes(list(set(self.idmap)), fields="symbol,name,taxid,go", species=self.species):

            if 'go' not in dic:
                continue
            if self.filter_go_category:
                dic['go'] = {category: dic['go'][category] for category in self.filter_go_category.split(",") if category in dic['go']}
            for category, go_list in dic['go'].items():
                if not isinstance(go_list, list):
                    go_list = [go_list]
                for go in go_list:
                    if (self.filter_go_evidence) and (go['evidence'] not in self.filter_go_evidence.split(",")):
                        continue
                    
                    current_info = {
                        'query': self.idmap[dic['_id']],
                        'mygene_id': dic['_id'],
                        'go_id': go['id'],
                        'go_term': go['term'],
                        'go_evidence': go['evidence'],
                        'go_category': category,
                        'symbol': dic['symbol'],
                        'name': dic['name'],
                        'taxid': dic['taxid']
                    }
                    info[self.idmap[dic['_id']]] = current_info
        return info
    
    def parse_go_based_info(self):
        """
        This converts the gene based info dictionary of dictionaries to a GO based info dictionary of lists.
        {go_id: [go_term, gene1, gene2, ...]}
        """
        # get unique go ids with the corresponding go terms
        info = {}
        for gene,dic in self.gene_based_info.items():
            info[dic['go_id']] = [dic['go_term']]

        # add all the genes associated to each go entry
        for gene,dic in self.gene_based_info.items():
            go_id = dic['go_id']
            if gene not in info[go_id]:
                info[go_id].append(gene)
        return info
    
    def save_to_tsv(self, filename: str) -> None:
        """
        It saves the parsed gene centric information in a tsv file.
        The parsing is performed in an sorted way following the query gene list order.
        """
        with open(filename, 'w') as f:
            f.write("\\t".join(self.gene_based_info[self.query[0]].keys()) + "\\n")
            for gene in self.query:  # sorted by query gene list
                if gene in self.gene_based_info:
                    f.write("\\t".join([str(val) for val in self.gene_based_info[gene].values()]) + "\\n")
    
    def save_to_gmt(self, filename: str) -> list:
        """
        It saves the parsed go centric information to a gmt file.
        The parsing is performed in an sorted way following the go id order.
        """
        info = dict(sorted(self.go_based_info.items(), key=lambda x: x[0]))
        with open(filename, 'w') as f:
            for go_id, go_list in info.items():
                tmp = sorted(go_list[1:])
                f.write(go_id + "\\t" + go_list[0] + "\\t" + "\\t".join(tmp) + "\\n")


if __name__ == "__main__":
    args = Arguments()

    # parse info for each gene
    mg = MyGene(args.input, args.species, args.filter_go_category, args.filter_go_evidence)

    # save info
    mg.save_to_tsv(args.output_tsv)
    mg.save_to_gmt(args.output_gmt)
    
    # write versions to file
    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, mygene])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))
