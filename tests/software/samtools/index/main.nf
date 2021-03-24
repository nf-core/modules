#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../../../../software/samtools/index/main.nf' addParams( options: [:] )

workflow test_samtools_index {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    SAMTOOLS_INDEX ( input )
}
