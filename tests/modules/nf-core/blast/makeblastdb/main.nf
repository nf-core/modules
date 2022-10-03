#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../../../../modules/blast/makeblastdb/main.nf'

workflow test_blast_makeblastdb {
    input =  [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
}
