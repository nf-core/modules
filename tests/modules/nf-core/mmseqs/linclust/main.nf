#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR          } from '../../../../../modules/nf-core/untar/main.nf'
include { MMSEQS_LINCLUST } from '../../../../../modules/nf-core/mmseqs/linclust/main.nf'

workflow test_mmseqs_linclust {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['mmseqs_tar_gz'], checkIfExists: true)
    ]

    ch_db = UNTAR ( input ).untar

    MMSEQS_LINCLUST ( ch_db )
}
