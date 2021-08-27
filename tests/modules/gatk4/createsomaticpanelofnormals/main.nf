#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../../modules/gatk4/createsomaticpanelofnormals/main.nf' addParams( options: [:] )

workflow test_gatk4_createsomaticpanelofnormals {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    GATK4_CREATESOMATICPANELOFNORMALS ( input )
}
