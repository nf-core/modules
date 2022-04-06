#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_LEARNREADORIENTATIONMODEL } from '../../../../modules/gatk4/learnreadorientationmodel/main.nf'

workflow test_gatk4_learnreadorientationmodel {

    input = [ [ id:'test' ], // meta map
              [file(params.test_data['homo_sapiens']['illumina']['test_test2_paired_mutect2_calls_f1r2_tar_gz'], checkIfExists: true)] ]

    GATK4_LEARNREADORIENTATIONMODEL ( input )
}
