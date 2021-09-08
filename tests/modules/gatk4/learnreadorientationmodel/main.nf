#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_LEARNREADORIENTATIONMODEL } from '../../../../modules/gatk4/learnreadorientationmodel/main.nf' addParams( options: [:] )

workflow test_gatk4_learnreadorientationmodel {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    GATK4_LEARNREADORIENTATIONMODEL ( input )
}
