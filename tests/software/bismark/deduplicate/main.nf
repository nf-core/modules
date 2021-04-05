#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_DEDUPLICATE } from '../../../../software/bismark/deduplicate/main.nf' addParams( options: [:] )

workflow test_bismark_deduplicate {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_methylated_paired_end_bam'], checkIfExists: true) ]
            ]

    BISMARK_DEDUPLICATE ( input )
}
