#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_MERGE } from '../../../../software/samtools/merge/main.nf' addParams( options: [:] )

workflow test_samtools_merge {

    def input = []
    input = [ [ id: 'test' ], // meta map 
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_methylated_paired_end.sorted.bam", checkIfExists: true ),
                file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/bam/test_single_end.sorted.bam", checkIfExists: true) ] ]

    SAMTOOLS_MERGE ( input )
}
