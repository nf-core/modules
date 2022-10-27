#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { BISMARK_DEDUPLICATE } from "$moduleDir/modules/nf-core/bismark/deduplicate/main.nf"

workflow test_bismark_deduplicate {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_bam'], checkIfExists: true) ]
            ]

    BISMARK_DEDUPLICATE ( input )
}
