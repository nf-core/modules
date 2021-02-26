#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT } from '../../../software/cnvkit/main.nf' addParams( options: [ 'args': '--method wgs --output-reference reference.cnn' ] )

workflow test_cnvkit {
    
    def fasta = file("${launchDir}/tests/data/fasta/human/human_subseq_chr21.fasta", checkIfExists: true)
    def annotationfile = file("${launchDir}/tests/data/txt/refflat.txt", checkIfExists: true)
       
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/bam/test_tumour_278_sub_chr21.bam", checkIfExists: true),
                file("${launchDir}/tests/data/bam/test_normal_280_sub_chr21.bam", checkIfExists: true) ] ]

    CNVKIT ( input, fasta, annotationfile )
}
