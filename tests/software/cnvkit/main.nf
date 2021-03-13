#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT } from '../../../software/cnvkit/main.nf' addParams( options: [ 'args': '--output-reference reference.cnn' ] )

workflow test_cnvkit {
    
    tumourbam = [ file("${launchDir}/tests/data/bam/test_paired_end.sorted_sarscov2.bam", checkIfExists: true) ]
    normalbam = [ file("${launchDir}/tests/data/bam/test_single_end.sorted_sarscov2.bam", checkIfExists: true) ]   
    
    def input = []
    input = [ [ id:'test' ], // meta map
              tumourbam, normalbam ]
    
    fasta = [ file("${launchDir}/tests/data/fasta/sarscov2/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true) ]
    targetfile = [ file("${launchDir}/tests/data/bed/baits.bed", checkIfExists: true) ]

    CNVKIT ( input, fasta, targetfile )
}
