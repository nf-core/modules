#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RTGTOOLS_VCFEVAL } from '../../../../modules/rtgtools/vcfeval/main.nf'

workflow test_rtgtools_vcfeval {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    RTGTOOLS_VCFEVAL ( input )
}
