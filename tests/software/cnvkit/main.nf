#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT } from '../../../software/cnvkit/main.nf' addParams( options: [ 'args': '--output-reference reference.cnn' ] )

workflow test_cnvkit {
    tumourbam = file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true)
    normalbam = file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_single_end.sorted.bam", checkIfExists: true)
    
    input = [ [ id:'test' ], // meta map
              tumourbam, 
              normalbam 
            ]
    fasta   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    targets = file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/baits.bed", checkIfExists: true)

    CNVKIT ( input, fasta, targets )
}
