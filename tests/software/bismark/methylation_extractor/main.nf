#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOME_PREPARATION    } from '../../../../software/bismark/genome_preparation/main.nf' addParams( options: [:] )
include { BISMARK_METHYLATION_EXTRACTOR } from '../../../../software/bismark/methylation_extractor/main.nf' addParams( options: [:] )

workflow test_bismark_methylation_extractor {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_methylated_paired_end.bam", checkIfExists: true) ] ]

    BISMARK_GENOME_PREPARATION ( file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) )

    BISMARK_METHYLATION_EXTRACTOR (
        input,
        BISMARK_GENOME_PREPARATION.out.index
    )
}
