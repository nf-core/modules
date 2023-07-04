#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR as UNTAR_QUERY   } from '../../../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR_TARGET  } from '../../../../../modules/nf-core/untar/main.nf'
include { MMSEQS_SEARCH          } from '../../../../../modules/nf-core/mmseqs/search/main.nf'

workflow test_mmseqs_search {


    input_query = [
        [ id:'test_query', single_end:true ], // meta map
        file(params.test_data['sarscov2']['genome']['mmseqs_tar_gz'], checkIfExists: true)
    ]
    ch_query_db = UNTAR_QUERY ( input_query  ).untar

    input_target = [
        [ id:'test_target', single_end:true ], // meta map
        file(params.test_data['sarscov2']['genome']['mmseqs_tar_gz'], checkIfExists: true)
    ]
    ch_target_db = UNTAR_TARGET ( input_target  ).untar

    MMSEQS_SEARCH ( ch_query_db, ch_target_db )
}
