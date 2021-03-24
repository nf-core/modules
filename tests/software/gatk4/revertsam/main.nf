#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_REVERTSAM } from '../../../../software/gatk4/revertsam/main.nf' addParams( options: [:] )

workflow test_gatk4_revertsam {
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.bam", checkIfExists: true) 
            ]

    GATK4_REVERTSAM ( input )
}
