#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTWGSMETRICS } from '../../../../software/picard/collectwgsmetrics/main.nf' addParams( options: [:] )

workflow test_picard_collectwgsmetrics {
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam.bai", checkIfExists: true)
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true) 

    PICARD_COLLECTWGSMETRICS ( input, fasta )
}
