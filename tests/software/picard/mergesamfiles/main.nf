#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_MERGESAMFILES } from '../../../../software/picard/mergesamfiles/main.nf' addParams( options: [:] )

workflow test_picard_mergesamfiles {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_single_end.bam", checkIfExists: true), ] 
            ]

    PICARD_MERGESAMFILES ( input )
}
