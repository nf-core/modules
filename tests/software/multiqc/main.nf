#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC } from '../../../software/multiqc/main.nf' addParams(options: [publish_dir: 'test_multiqc'])
include { test_paired_end } from '../fastqc/main.nf'  addParams( options: [ publish_dir:'test_paired_end' ] )

workflow test_multiqc {
    test_paired_end()


//    test_paired_end.out.html.collect { it[1] }
//        .join(test_paired_end.out.zip)
//        .view()

    input = [
        test_paired_end.out.zip.collect { it[1] }.ifEmpty([])
    ]

    MULTIQC(input)
}
