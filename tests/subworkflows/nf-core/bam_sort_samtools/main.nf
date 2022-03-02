#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_SORT_SAMTOOLS } from '../../../../subworkflows/nf-core/bam_sort_samtools/main' addParams( sort_options: ['suffix': '.sorted']  )

workflow test_bam_sort_samtools_single_end {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
            ]

     BAM_SORT_SAMTOOLS ( input )
}

workflow test_bam_sort_samtools_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]

     BAM_SORT_SAMTOOLS ( input )
}
