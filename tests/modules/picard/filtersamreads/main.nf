#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_FILTERSAMREADS } from '../../../../modules/picard/filtersamreads/main.nf' addParams( options: [:] )

workflow test_picard_filtersamreads {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
    filter = 'includeAligned'

    PICARD_FILTERSAMREADS ( input, filter, [] )
}

workflow test_picard_filtersamreads_readlist {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
    filter = 'includeReadAligned'
    #readlist = file(params.test_data['sarscov2']['illumina']['picard_readlist_test_single_end_bam'], checkIfExists: true) ]
    readlist = file('https://github.com/jfy133/test-datasets/raw/modules/data/genomics/sarscov2/illumina/bam/test.single_end.bam.alignedreadlist.txt', checkIfExists: true)

    PICARD_FILTERSAMREADS ( input, filter, readlist )
}
