#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR      } from '../../../../software/untar/main.nf'      addParams( options: [:] )
include { LAST_TRAIN } from '../../../../software/last/train/main.nf' addParams( options: [:] )

workflow test_last_train {

    db =         [ file(params.test_data['sarscov2']['genome']['lastdb_tar_gz'], checkIfExists: true) ]
    input =      [ [ id:'contigs' ], // meta map
                   file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    UNTAR ( db )
    // add meta map needed by LAST_TRAIN
    ch_index = UNTAR.out.untar.map {row -> [ [ id : 'genome' ], row ] }
    LAST_TRAIN ( ch_index, input )
}
