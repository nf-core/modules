#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOMEPREPARATION    } from '../../../../software/bismark/genomepreparation/main.nf' addParams( options: [:] )
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../../software/bismark/methylationextractor/main.nf' addParams( options: [:] )

workflow test_bismark_methylationextractor {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_methylated_paired_end.bam", checkIfExists: true) ] ]

    BISMARK_GENOMEPREPARATION ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )

    BISMARK_METHYLATIONEXTRACTOR (
        input,
        BISMARK_GENOMEPREPARATION.out.index
    )
}
