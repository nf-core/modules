#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC } from '../../../software/multiqc/main.nf' addParams(options: [publish_dir: 'test_multiqc'])
include { test_paired_end } from '../fastqc/main.nf'  addParams( options: [ publish_dir:'test_paired_end' ] )

workflow test_multiqc {
    test_paired_end()

    input = [
        [id: 'test'],
        test_paired_end.out.html.join(test_paired_end.out.zip, by: 0).flatten().filter(java.nio.file.Path).toList()
    ]

    MULTIQC(input)
}
