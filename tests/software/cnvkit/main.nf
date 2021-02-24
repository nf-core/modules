#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT } from '../../../software/cnvkit/main.nf' addParams( options: [ 'args': '--targets my_baits.bed --method wgs --annotate refflat.txt --output-reference reference.cnn' ] )

workflow test_cnvkit {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/bam/test_tumour_278_sub_chr21.bam", checkIfExists: true),
                file("${launchDir}/tests/data/bam/test_normal_280_sub_chr21.bam", checkIfExists: true),
                file("${launchDir}/tests/data/fasta/human/human_subseq_chr21.fasta", checkIfExists: true) ] ]

    CNVKIT ( input )
}
