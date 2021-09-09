#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_LEARNREADORIENTATIONMODEL } from '../../../../modules/gatk4/learnreadorientationmodel/main.nf' addParams( options: [:] )

workflow test_gatk4_learnreadorientationmodel {

    input = [ [ id:'test', single_end:false ], // meta map
              [file("/home/AD/gmackenz/test_data/test.f1r2.tar.gz", checkIfExists: true)] ]

    GATK4_LEARNREADORIENTATIONMODEL ( input )
}
