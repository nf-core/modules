#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MARKDUPLICATES } from '../../../../software/gatk4/markduplicates/main.nf' addParams( options: [:] )

workflow test_gatk4_markduplicates {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    GATK4_MARKDUPLICATES ( input )
}
