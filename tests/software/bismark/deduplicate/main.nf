#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_DEDUPLICATE } from '../../../../software/bismark/deduplicate/main.nf' addParams( options: [:] )

workflow test_bismark_deduplicate {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe.bam", checkIfExists: true) ] ]

    BISMARK_DEDUPLICATE ( input )
}
