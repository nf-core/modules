#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MMSEQS_SEARCH } from '../../../../../modules/nf-core/mmseqs/search/main.nf'

workflow test_mmseqs_search {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    MMSEQS_SEARCH ( input )
}
