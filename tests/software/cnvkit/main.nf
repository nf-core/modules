#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT } from '../../../../software/cnvkit/main.nf' addParams( options: [ publish_dir:'test_with_normal' ] )
include { CNVKIT } from '../../../../software/cnvkit/main.nf' addParams( options: [ publish_dir:'test_without_normal' ] )

/*
 * Test with normal sample
 */
workflow test_with_normal {

    def input = []
    input = [ [ id:'test', with_normal:true ], // meta map
              [ file("${launchDir}/tests/data/bam/tumour_278_sub_chr21.bam", checkIfExists: true),
	        file("${launchDir}/tests/data/bam/normal_280_sub_chr21.bam", checkIfExists: true) ] ]

    CNVKIT (
        input,
        file("${launchDir}/tests/data/fasta/human/human_subseq_chr21.fasta", checkIfExists: true),
        file("${launchDir}/tests/data/annotation/refflat.txt", checkIfExists: true)
    )
}

/*
 * Test without normal sample
 */
workflow test_without_normal {

    def input = []
    input = [ [ id:'test', with_normal:false ], // meta map
              [ file("${launchDir}/tests/data/bam/tumour_278_sub_chr21.bam", checkIfExists: true) ] ]

    CNVKIT (
        input,
        file("${launchDir}/tests/data/fasta/human/human_subseq_chr21.fasta", checkIfExists: true),
        file("${launchDir}/tests/data/annotation/refflat.txt", checkIfExists: true)
    )
}

