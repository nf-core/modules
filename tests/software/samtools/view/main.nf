#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_VIEW } from '../../../../software/samtools/view/main.nf' addParams( options: [:] )

workflow test_samtools_view {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)

            ]

    SAMTOOLS_VIEW ( input )
}
