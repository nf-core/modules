#!/usr/bin/env nextflow



include { BISMARK_DEDUPLICATE } from '../../../../modules/bismark/deduplicate/main.nf'

workflow test_bismark_deduplicate {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_bam'], checkIfExists: true) ]
            ]

    BISMARK_DEDUPLICATE ( input )
}
