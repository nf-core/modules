# TO be imported from refget package when it is finished
# from refget import fasta_checksum

import binascii
import hashlib
import os

import pyfaidx


def trunc512_digest(seq, offset=24):
    digest = hashlib.sha512(seq.encode()).digest()
    hex_digest = binascii.hexlify(digest[:offset])
    return str(hex_digest.decode())


def parse_fasta(fa_file):
    try:
        fa_object = pyfaidx.Fasta(fa_file)
    except pyfaidx.UnsupportedCompressionFormat:
        # pyfaidx can handle bgzip but not gzip; so we just hack it here and
        # unzip the file for checksumming, then rezip it for the rest of the
        # asset build.
        # TODO: streamline this to avoid repeated compress/decompress
        # in refgenie we feed this function with uncompressed, newly built
        # FASTA file, so compression issues are not relevant
        os.system("gunzip {}".format(fa_file))
        fa_file_unzipped = fa_file.replace(".gz", "")
        fa_object = pyfaidx.Fasta(fa_file_unzipped)
        os.system("gzip {}".format(fa_file_unzipped))
    return fa_object


def fasta_checksum(fa_file, checksum_function=trunc512_digest):
    """
    Just calculate checksum of fasta file without loading it.
    """
    fa_object = parse_fasta(fa_file)
    content_checksums = {}
    for k in fa_object.keys():
        content_checksums[k] = checksum_function(str(fa_object[k]))
    collection_string = ";".join([":".join(i) for i in content_checksums.items()])
    collection_checksum = checksum_function(collection_string)
    return collection_checksum, content_checksums
