#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { UNTAR as UNTAR_QUERY   } from '../../../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR_QUERY2  } from '../../../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR_TARGET  } from '../../../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR_RESULT  } from '../../../../../modules/nf-core/untar/main.nf'
include { UNTAR as UNTAR_RESULT2 } from '../../../../../modules/nf-core/untar/main.nf'
include { MMSEQS_CREATETSV       } from '../../../../../modules/nf-core/mmseqs/createtsv/main.nf'

workflow test_mmseqs_createtsv {

    input_file = file(params.test_data['sarscov2']['genome']['mmseqs_tar_gz'], checkIfExists: true)

    ch_query = UNTAR_QUERY ([[ id:'test_query', single_end:true ],input_file] ).untar
    ch_target = UNTAR_TARGET( [[ id:'test_target', single_end:true ],input_file] ).untar
    ch_result = UNTAR_RESULT( [[ id:'test_result', single_end:true ],input_file] ).untar

    MMSEQS_CREATETSV (ch_result, ch_query, ch_target  )
}

workflow test_mmseqs_createtsv_optional {

    input_file = file(params.test_data['sarscov2']['genome']['mmseqs_tar_gz'], checkIfExists: true)

    ch_query = UNTAR_QUERY2 ([[ id:'test_query', single_end:true ],input_file] ).untar
    ch_result = UNTAR_RESULT2( [[ id:'test_result', single_end:true ],input_file] ).untar

    MMSEQS_CREATETSV (ch_result, ch_query, [[:],[]]  )
}
