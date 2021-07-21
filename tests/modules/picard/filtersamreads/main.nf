#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_SORTSAM        } from '../../../../modules/picard/sortsam/main.nf'        addParams( options: [:] )
include { PICARD_FILTERSAMREADS } from '../../../../modules/picard/filtersamreads/main.nf' addParams( options: [:] )

workflow test_picard_filtersamreads {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
    sort_order = 'queryname'
    filter = 'includeAligned'

    PICARD_SORTSAM ( input, sort_order  )
    PICARD_FILTERSAMREADS ( PICARD_SORTSAM.out.bam, filter, [] )
}

workflow test_picard_filtersamreads_readlist {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
    filter = 'includeReadList'
    readlist = file(params.test_data['sarscov2']['illumina']['picard_readlist_test_single_end_bam'], checkIfExists: true)

    PICARD_FILTERSAMREADS ( input, filter, readlist )
}
