#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR      } from '../../../../modules/untar/main.nf'      addParams( options: [:] )
include { LAST_TRAIN } from '../../../../modules/last/train/main.nf' addParams( options: [:] )

workflow test_last_train {

    db =         [ file(params.test_data['sarscov2']['genome']['lastdb_tar_gz'], checkIfExists: true) ]
    input =      [ [ id:'contigs' ], // meta map
                   file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true) ]
    UNTAR ( db )
    LAST_TRAIN ( input, UNTAR.out.untar )
}
