#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_FILTERSAMREADS } from '../../../../modules/picard/filtersamreads/main.nf' addParams( options: [:] )

workflow test_picard_filtersamreads {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    PICARD_FILTERSAMREADS ( input )
}
