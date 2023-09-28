#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../../../../../modules/nf-core/blast/makeblastdb/main.nf'
include { BLAST_BLASTP } from '../../../../../modules/nf-core/blast/blastp/main.nf'

workflow test_blast_blastp {

    input = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
    BLAST_BLASTP ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db )
}
