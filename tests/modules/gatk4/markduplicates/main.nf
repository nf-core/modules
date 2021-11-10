#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MARKDUPLICATES } from '../../../../modules/gatk4/markduplicates/main.nf' addParams( options: [:] )

workflow test_gatk4_markduplicates_keep_metrics {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    GATK4_MARKDUPLICATES ( input, true )
}

workflow test_gatk4_markduplicates_no_metrics {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    GATK4_MARKDUPLICATES ( input, false )
}

workflow test_gatk4_markduplicates_multiple_bams_keep_metrics {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
            ]

    GATK4_MARKDUPLICATES ( input, true )
}
