import binascii
import hashlib
import logging
import os
from gzip import open as gzopen

from .exceptions import RefgenconfError
from .henge import ITEM_TYPE, Henge


def trunc512_digest(seq, offset=24):
    digest = hashlib.sha512(seq.encode()).digest()
    hex_digest = binascii.hexlify(digest[:offset])
    return hex_digest.decode()


# module constants
def _schema_path(name):
    return os.path.join(SCHEMA_FILEPATH, name)


CONTENT_ALL_A_IN_B = 2 ** 0
CONTENT_ALL_B_IN_A = 2 ** 1
LENGTHS_ALL_A_IN_B = 2 ** 2
LENGTHS_ALL_B_IN_A = 2 ** 3
NAMES_ALL_A_IN_B = 2 ** 4
NAMES_ALL_B_IN_A = 2 ** 5
CONTENT_A_ORDER = 2 ** 6
CONTENT_B_ORDER = 2 ** 7
CONTENT_ANY_SHARED = 2 ** 8
LENGTHS_ANY_SHARED = 2 ** 9
NAMES_ANY_SHARED = 2 ** 10

FLAGS = {
    CONTENT_ALL_A_IN_B: "CONTENT_ALL_A_IN_B",
    CONTENT_ALL_B_IN_A: "CONTENT_ALL_B_IN_A",
    LENGTHS_ALL_A_IN_B: "LENGTHS_ALL_A_IN_B",
    LENGTHS_ALL_B_IN_A: "LENGTHS_ALL_B_IN_A",
    NAMES_ALL_A_IN_B: "NAMES_ALL_A_IN_B",
    NAMES_ALL_B_IN_A: "NAMES_ALL_B_IN_A",
    CONTENT_ANY_SHARED: "CONTENT_ANY_SHARED",
    LENGTHS_ANY_SHARED: "LENGTHS_ANY_SHARED",
    NAMES_ANY_SHARED: "NAMES_ANY_SHARED",
    CONTENT_A_ORDER: "CONTENT_A_ORDER",
    CONTENT_B_ORDER: "CONTENT_B_ORDER",
}

NAME_KEY = "name"
SEQ_KEY = "sequence"
LEN_KEY = "length"

# internal schemas paths determination
ASL_NAME = "AnnotatedSequenceList"
ASDL_NAME = "AnnotatedSequenceDigestList"
SCHEMA_NAMES = [ASL_NAME, ASDL_NAME]
SCHEMA_FILEPATH = os.path.join(os.path.dirname(__file__), "schemas")
INTERNAL_SCHEMAS = [_schema_path(f"{s}.yaml") for s in SCHEMA_NAMES]

_LOGGER = logging.getLogger(__name__)


class SeqColClient(Henge):
    """
    Extension of henge that accommodates collections of sequences.
    """

    def __init__(
        self, database, schemas=None, henges=None, checksum_function=trunc512_digest
    ):
        """
        A user interface to insert and retrieve decomposable recursive unique
        identifiers (DRUIDs).

        :param dict database: Dict-like lookup database with sequences
            and hashes
        :param dict schemas: One or more jsonschema schemas describing the
            data types stored by this Henge
        :param function(str) -> str checksum_function: Default function to
            handle the digest of the
            serialized items stored in this henge.
        """
        assert all([os.path.exists(s) for s in INTERNAL_SCHEMAS]), RefgenconfError(
            f"Missing schema files: {INTERNAL_SCHEMAS}"
        )
        super(SeqColClient, self).__init__(
            database=database,
            schemas=schemas or INTERNAL_SCHEMAS,
            henges=henges,
            checksum_function=checksum_function,
        )

    def load_fasta(self, fa_file, skip_seq=False, gzipped=False):
        """
        Load a sequence collection into the database

        :param str fa_file: path to the FASTA file to parse and load
        :param bool skip_seq: whether to disregard the actual sequences,
            load just the names and lengths
        :param bool skip_seq: whether to disregard the actual sequences,
            load just the names and lengths
        """
        seq = ""
        name = ""
        init = False
        aslist = []
        openfun = gzopen if gzipped else open
        with openfun(fa_file, "rt") as f:
            for line in f:
                line = line.strip("\n")
                if line.startswith(">"):
                    if not init:
                        name = line.replace(">", "")
                    else:
                        aslist.append(
                            {
                                NAME_KEY: name,
                                LEN_KEY: len(seq),
                                SEQ_KEY: "" if skip_seq else trunc512_digest(seq),
                            }
                        )
                        name = line.replace(">", "")
                    seq = ""
                    continue
                init = True
                seq = seq + line
            aslist.append(
                {
                    NAME_KEY: name,
                    LEN_KEY: len(seq),
                    SEQ_KEY: "" if skip_seq else trunc512_digest(seq),
                }
            )

        collection_checksum = self.insert(aslist, ASDL_NAME)
        _LOGGER.info(f"Loaded {ASDL_NAME} ({len(aslist)} sequences)")
        return collection_checksum, aslist

    @staticmethod
    def compare_asds(asdA, asdB, explain=False):
        """
        Compare Annotated Sequence Digests (ASDs) -- digested sequences and metadata

        :param str asdA: ASD for first sequence collection to compare.
        :param str asdB: ASD for second sequence collection to compare.
        :param bool explain: Print an explanation of the flag? [Default: False]
        """

        def _xp(prop, lst):
            """Extract property from a list of dicts"""
            return list(map(lambda x: x[prop], lst))

        def _index(x, lst):
            """Find an index of a sequence element in a list of dicts"""
            try:
                return _xp(SEQ_KEY, lst).index(x)
            except:
                return None

        def _get_common_content(lstA, lstB):
            """
            Find the intersection between two list of dicts with sequences
            """
            return list(
                filter(None.__ne__, [_index(x, lstB) for x in _xp(SEQ_KEY, lstA)])
            )

        # Not ideal, but we expect these to return lists, but if the item was
        # singular only a dict is returned
        if not isinstance(asdA, list):
            asdA = [asdA]
        if not isinstance(asdB, list):
            asdB = [asdB]

        ainb = [x in _xp(SEQ_KEY, asdB) for x in _xp(SEQ_KEY, asdA)]
        bina = [x in _xp(SEQ_KEY, asdA) for x in _xp(SEQ_KEY, asdB)]

        return_flag = 0  # initialize
        if any(ainb):
            ordA = _get_common_content(asdA, asdB)
            if ordA == sorted(ordA):
                return_flag += CONTENT_A_ORDER
        if any(bina):
            ordB = _get_common_content(asdB, asdA)
            if ordB == sorted(ordB):
                return_flag += CONTENT_B_ORDER

        ainb_len = [x in _xp(LEN_KEY, asdB) for x in _xp(LEN_KEY, asdA)]
        bina_len = [x in _xp(LEN_KEY, asdA) for x in _xp(LEN_KEY, asdB)]

        ainb_name = [x in _xp(NAME_KEY, asdB) for x in _xp(NAME_KEY, asdA)]
        bina_name = [x in _xp(NAME_KEY, asdA) for x in _xp(NAME_KEY, asdB)]

        if any(ainb):
            return_flag += CONTENT_ANY_SHARED
        if all(ainb):
            return_flag += CONTENT_ALL_A_IN_B
        if all(bina):
            return_flag += CONTENT_ALL_B_IN_A

        if any(ainb_name):
            return_flag += NAMES_ANY_SHARED
        if all(ainb_name):
            return_flag += NAMES_ALL_A_IN_B
        if all(bina_name):
            return_flag += NAMES_ALL_B_IN_A

        if any(ainb_len):
            return_flag += LENGTHS_ANY_SHARED
        if all(ainb_len):
            return_flag += LENGTHS_ALL_A_IN_B
        if all(bina_len):
            return_flag += LENGTHS_ALL_B_IN_A

        if explain:
            explain_flag(return_flag)
        return return_flag

    def compare(self, digestA, digestB, explain=False):
        """
        Given two collection checksums in the database, provide some information
        about how they are related.

        :param str digestA: Digest for first sequence collection to compare.
        :param str digestB: Digest for second sequence collection to compare.
        :param bool explain: Print an explanation of the flag? [Default: False]
        """
        typeA = self.database[digestA + ITEM_TYPE]
        typeB = self.database[digestB + ITEM_TYPE]

        if typeA != typeB:
            _LOGGER.error(
                f"Can't compare objects of different types: " f"{typeA} vs {typeB}"
            )

        asdA = self.retrieve(digestA, reclimit=1)
        asdB = self.retrieve(digestB, reclimit=1)
        return self.compare_asds(asdA, asdB, explain=explain)


# Static functions below (these don't require a database)


def explain_flag(flag):
    """Explains a compare flag"""
    print(f"Flag: {flag}\nBinary: {bin(flag)}\n")
    for e in range(0, 13):
        if flag & 2 ** e:
            print(FLAGS[2 ** e])
