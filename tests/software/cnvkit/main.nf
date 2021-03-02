#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT } from '../../../software/cnvkit/main.nf' addParams( options: [ 'args': '--method wgs --output-reference reference.cnn' ] )

workflow test_cnvkit {
    
    tumourbam = [ file("${launchDir}/tests/data/bam/test_tumour_278_sub_chr21.bam", checkIfExists: true) ]
    normalbam = [ file("${launchDir}/tests/data/bam/test_normal_280_sub_chr21.bam", checkIfExists: true) ]   
    
    def input = []
    input = [ [ id:'test' ], // meta map
              tumourbam, normalbam ]
    
    fasta = [ file("${launchDir}/tests/data/fasta/human/human_subseq_chr21.fasta", checkIfExists: true) ]
    annotationfile = [ file("${launchDir}/tests/data/txt/refflat.txt", checkIfExists: true) ]

    CNVKIT ( input, fasta, annotationfile )
}
