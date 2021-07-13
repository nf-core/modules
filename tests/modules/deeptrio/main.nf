#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPTRIO } from '../../../modules/deeptrio/main.nf' addParams( options: [:] )

workflow test_deeptrio {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    DEEPTRIO ( input )
}
