#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT } from '../../../software/cnvkit/main.nf' addParams( options: [ 'args': '--output-reference reference.cnn' ] )

workflow test_cnvkit {
    
    tumourbam = [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ]
    normalbam = [ file("${launchDir}/tests/data/genomics/sarscov2/bam/test_single_end.sorted.bam", checkIfExists: true) ]   
    
    def input = []
    input = [ [ id:'test' ], // meta map
              tumourbam, normalbam ]
    
    fasta = [ file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ]
    targetfile = [ file("${launchDir}/tests/data/genomics/sarscov2/bed/baits.bed", checkIfExists: true) ]

    CNVKIT ( input, fasta, targetfile )
}
