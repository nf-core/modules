#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_SORTSAM        } from '../../../../modules/picard/sortsam/main.nf'
include { PICARD_FILTERSAMREADS } from '../../../../modules/picard/filtersamreads/main.nf'

workflow test_picard_filtersamreads {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
    sort_order = 'queryname'
    filter = 'includeAligned'

    PICARD_SORTSAM ( input, sort_order )
    PICARD_SORTSAM.out.bam
        .map {
            [ it[0], it[1], [] ]
        }
        .set{ ch_sorted_for_filtersamreads }
    PICARD_FILTERSAMREADS ( ch_sorted_for_filtersamreads, filter )
}

workflow test_picard_filtersamreads_readlist {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam_readlist_txt'], checkIfExists: true) ]
    filter = 'includeReadList'

    PICARD_FILTERSAMREADS ( input, filter )
}
