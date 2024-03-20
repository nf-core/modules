#!/usr/bin/env python3
import argparse
import mygene
import shlex


"""
This python script uses the mygene module to query the MyGene.info API and
retrieve the go terms associated with a list of gene ids. The gene ids should
ideally be Ensembl or Entrez ids. The script generates two outputs:
    1. A tsv file containing information related to each query. The columns
        include query, mygene_id, go_id, go_term, go_evidence, go_category,
        symbol, name, and taxid.
    2. A gmt file containing information related to each go term. Each row
        includes the go id, go term, and the genes associated with that go term.

Author: Suzanne Jin
License: Apache 2.0 (same as the mygene library)
"""


class Arguments:
    """
    Parses the argments, including the ones coming from $task.ext.args.
    """
    def __init__(self) -> None:
        self.input = "$gene_list"
        self.prefix = (
            "$task.ext.prefix"
            if "$task.ext.prefix" != "null"
            else "$meta.id"
        )
        self.output_gmt = self.prefix + ".gmt"
        self.output_tsv = self.prefix + ".tsv"
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
        # input parameters
        parser.add_argument('--columname', default='gene_id', help='Name of the column where the gene ids are stored in the input file. Default: gene_id')
        # filtering parameters
        parser.add_argument('--species', default=None, help="Comma separated of common name of the species or taxon ids")
        parser.add_argument('--go_category', default=None, help="Comma separated list of GO categories to keep. Default: all")
        parser.add_argument('--go_evidence', default=None, help="Comma separated list of GO evidence codes to keep. Default: all")
        # additional parameters for querymany
        parser.add_argument('--scopes', default=None, help="Comma separated list of scopes to search for.")
        parser.add_argument('--entrezonly', default=False, help="When true, the query returns only the hits with valid Entrez gene ids. Default: false.")
        parser.add_argument('--ensemblonly', default=False, help="When true, the query returns only the hits with valid Ensembl gene ids. Default: False")
        # output parameters
        parser.add_argument('--generate_tsv', default=False, help="Also generate a tsv file with the gene based information. Default: False")
        args = parser.parse_args(args_list)

        # Convert "null" values to default values
        # convert true to True and false to False
        for attr in vars(args):
            value = getattr(args, attr)
            if value == "null":
                setattr(args, attr, parser.get_default(attr))
            elif value == "true":
                setattr(args, attr, True)
            elif value == "false":
                setattr(args, attr, False)

        # check if the arguments are valid
        if args.go_category:
            args.go_category = args.go_category.upper()
            for category in args.go_category.split(","):
                if category not in ["BP", "MF", "CC"]:
                    raise ValueError("The GO category should be one of BP, MF, or CC.")
        if args.go_evidence:
            args.go_evidence = args.go_evidence.upper()

        # Assign args attributes to self attributes
        for attr in vars(args):
            setattr(self, attr, getattr(args, attr))

    def print_args(self) -> None:
        """
        Print the arguments.
        """
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")


class Version:
    """
    Parse the versions of the modules used in the script.
    """

    @staticmethod
    def get_versions(modules: list) -> dict:
        """
        This function takes a list of modules and returns a dictionary with the
        versions of each module.
        """
        return {module.__name__: module.__version__ for module in modules}

    @staticmethod
    def format_yaml_like(data: dict, indent: int = 0) -> str:
        """
        Formats a dictionary to a YAML-like string.

        Args:
            data (dict): The dictionary to format.
            indent (int): The current indentation level.

        Returns:
            yaml_str: A string formatted as YAML.
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
    This class will query the MyGene.info API and retrieve the go terms
    associated with a list of gene ids.

    Args:
        query (list) :
            A list of gene ids. Ideally Ensembl or Entrez ids.
        species (str) :
            Comma separated of common name of the species or taxon ids.
            If not provided, the API will return information for all species.
        go_category (str) :
            Comma separated list of GO categories to keep. If not provided,
            the API will return all categories.
        go_evidence (str) :
            Comma separated list of GO evidence codes to keep. If not provided,
            the API will return all evidence codes.
    """
    def __init__(self, query: list, species: str, scopes: str, entrezonly: bool, ensemblonly: bool, go_category: str = None, go_evidence: str = None) -> None:
        self.query = query
        self.fields = "go,symbol,name,taxid"
        self.species = species
        self.scopes = scopes
        self.entrezonly = entrezonly
        self.ensemblonly = ensemblonly
        self.go_category = go_category
        self.go_evidence = go_evidence
        self.mg = mygene.MyGeneInfo()

        # query info
        self.idmap = self.query2idmap()
        print(f"fetched {len(self.idmap)} ids from {len(self.query)} queries")

        # fetch and parse go information into the correct format
        self.gene_based_info = self.parse_gene_based_info()
        self.go_based_info = self.parse_go_based_info()
        print(f"parsed {len(self.gene_based_info)} gene centric info and {len(self.go_based_info)} go centric info")

    def query2idmap(self) -> dict:
        """
        It returns a dictionary with the query ids as keys and the mygene ids as values.
        """
        q = self.mg.querymany(
            self.query,
            scopes=self.scopes,
            species=self.species,
            entrezonly=self.entrezonly,
            ensemblonly=self.ensemblonly,
            returnall=True
        )
        return {dic['_id']: dic['query'] for dic in q['out'] if '_id' in dic}

    def id2info(self) -> list:
        """
        It returns a list of dictionaries with the info returned from getgenes for all the query ids.
        """
        return self.mg.getgenes(list(set(self.idmap)), fields=self.fields, species=self.species)

    def parse_gene_based_info(self) -> dict:
        """
        It parses the information for all the query ids.
        At the end it returns a dictionary {query gene: {}} of dictionaries with the following keys:
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
        for dic in self.id2info():

            if 'go' not in dic:
                continue
            if self.go_category:
                dic['go'] = {category: dic['go'][category] for category in self.go_category.split(",") if category in dic['go']}
            for category, go_list in dic['go'].items():
                if not isinstance(go_list, list):
                    go_list = [go_list]
                for go in go_list:
                    if (self.go_evidence) and (go['evidence'] not in self.go_evidence.split(",")):
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
            if gene in info[go_id]:
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


def load_list(filename: str, columname: str) -> list:
        """
        It loads the list of gene ids from a file.
        The columname is the name of the column where the gene ids are stored.
        """
        if filename.split('.')[-1] == 'tsv':
            sep = "\\t"
        elif filename.split('.')[-1] == 'csv':
            sep = ","
        else:
            raise ValueError("The input file extension should be either tsv or csv.")
        with open(filename, 'r') as f:
            idx = f.readline().strip().split(sep).index(columname)
            return [line.strip().split(sep)[idx] for line in f]


if __name__ == "__main__":

    # parse and print arguments
    args = Arguments()
    args.print_args()

    # load gene list
    gene_list = load_list(args.input, args.columname)

    # run mygene api
    mg = MyGene(
            gene_list,
            species=args.species,
            scopes=args.scopes,
            entrezonly=args.entrezonly,
            ensemblonly=args.ensemblonly,
            go_category=args.go_category,
            go_evidence=args.go_evidence
        )

    # save info
    mg.save_to_gmt(args.output_gmt)
    if args.generate_tsv:
        mg.save_to_tsv(args.output_tsv)

    # write versions to file
    versions_this_module = {}
    versions_this_module["${task.process}"] = Version.get_versions([argparse, mygene])
    with open("versions.yml", "w") as f:
        f.write(Version.format_yaml_like(versions_this_module))
