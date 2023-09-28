#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../../../../../modules/nf-core/blast/makeblastdb/main.nf'
include { BLAST_TBLASTN } from '../../../../../modules/nf-core/blast/tblastn/main.nf'

workflow test_blast_tblastn {

    input = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    input_pep = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
    BLAST_TBLASTN ( [ [id:'test'], input_pep ], BLAST_MAKEBLASTDB.out.db )
}
