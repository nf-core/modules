#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_LEARNREADORIENTATIONMODEL } from '../../../../modules/gatk4/learnreadorientationmodel/main.nf' addParams( options: [:] )

workflow test_gatk4_learnreadorientationmodel {

    input = [ [ id:'test' ], // meta map
              [file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.f1r2.tar.gz", checkIfExists: true)] ]

    GATK4_LEARNREADORIENTATIONMODEL ( input )
}
