#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC } from '../../../software/multiqc/main.nf' addParams(options: [publish_dir: 'test_multiqc'])
include { test_paired_end } from '../fastqc/main.nf' addParams(options: [publish_dir: 'test_paired_end'])

workflow test_multiqc {
    test_paired_end()

    input = [
        test_paired_end.out.zip.collect { it[1] }.ifEmpty([])
    ]

    MULTIQC(*input)
}
