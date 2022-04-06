#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../modules/diamond/makedb/main.nf'
include { DIAMOND_BLASTX } from '../../../../modules/diamond/blastx/main.nf'

workflow test_diamond_blastx {

    db = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    fasta = [ file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]

    DIAMOND_MAKEDB ( db )
    DIAMOND_BLASTX ( [ [id:'test'], fasta ], DIAMOND_MAKEDB.out.db )
}
