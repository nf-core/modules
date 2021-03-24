#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_DEDUPLICATE } from '../../../../software/bismark/deduplicate/main.nf' addParams( options: [:] )

workflow test_bismark_deduplicate {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_methylated_paired_end.bam", checkIfExists: true) ]
            ]

    BISMARK_DEDUPLICATE ( input )
}
