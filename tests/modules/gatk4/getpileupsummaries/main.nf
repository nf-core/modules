#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GETPILEUPSUMMARIES } from '../../../../modules/gatk4/getpileupsummaries/main.nf' addParams( options: [:] )

workflow test_gatk4_getpileupsummaries {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    GATK4_GETPILEUPSUMMARIES ( input )
}
