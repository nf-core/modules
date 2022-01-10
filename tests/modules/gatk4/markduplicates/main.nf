#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MARKDUPLICATES } from '../../../../modules/gatk4/markduplicates/main.nf'

workflow test_gatk4_markduplicates {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    GATK4_MARKDUPLICATES ( input )
}

workflow test_gatk4_markduplicates_multiple_bams {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
            ] ]

    GATK4_MARKDUPLICATES ( input )
}
