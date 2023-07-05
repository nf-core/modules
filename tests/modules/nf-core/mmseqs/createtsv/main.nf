#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { MMSEQS_CREATETSV      } from '../../../../../modules/nf-core/mmseqs/createtsv/main.nf'
include { MMSEQS_CREATEDB      } from '../../../../../modules/nf-core/mmseqs/createtsv/main.nf'

workflow test_mmseqs_createtsv {

    input_file = file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)

    ch_query = MMSEQS_CREATEDB ([[ id:'test_query', single_end:true ],input_file] ).db
    ch_target = MMSEQS_CREATEDB( [[ id:'test_target', single_end:true ],input_file] ).untar
    ch_result = UNTAR_RESULT( [[ id:'test_result', single_end:true ],input_file] ).untar

    MMSEQS_CREATETSV (ch_result, ch_query, ch_target  )
}
