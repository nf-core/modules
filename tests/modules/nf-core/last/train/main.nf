#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR      } from '../../../../../modules/nf-core/untar/main.nf'
include { LAST_TRAIN } from '../../../../../modules/nf-core/last/train/main.nf'

workflow test_last_train {

    db =         [ [], file(params.test_data['sarscov2']['genome']['lastdb_tar_gz'], checkIfExists: true) ]
    input =      [ [ id:'contigs' ], // meta map
                   file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    UNTAR ( db )
    LAST_TRAIN ( input, UNTAR.out.untar.map{ it[1] } )
}
