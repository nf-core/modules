#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_MARKDUPLICATES } from '../../../../software/picard/markduplicates/main.nf' addParams( options: [:] )

workflow test_picard_markduplicates_sorted_bam  {
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) 
            ]

    PICARD_MARKDUPLICATES ( input )
}

workflow test_picard_markduplicates_unsorted_bam  {
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.bam", checkIfExists: true) 
            ]

    PICARD_MARKDUPLICATES ( input )
}
