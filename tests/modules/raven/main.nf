#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RAVEN } from '../../../modules/raven/main.nf' addParams( options: [:] )

workflow test_raven {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    RAVEN ( input )
}
