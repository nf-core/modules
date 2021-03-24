#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTWGSMETRICS } from '../../../../software/picard/collectwgsmetrics/main.nf' addParams( options: [:] )

workflow test_picard_collectwgsmetrics {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true),
            file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam.bai", checkIfExists: true),]

    PICARD_COLLECTWGSMETRICS ( 
        input,
        file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) 
        )
}
