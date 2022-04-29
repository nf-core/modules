#!/usr/bin/env nextflow



include { BLAST_MAKEBLASTDB } from '../../../../modules/blast/makeblastdb/main.nf'

workflow test_blast_makeblastdb {
    input =  [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
}
