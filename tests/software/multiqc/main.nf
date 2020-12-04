#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC } from '../../../software/multiqc/main.nf' addParams(options: [publish_dir: 'test_multiqc'])
include { test_paired_end } from '../fastqc/main.nf' addParams(options: [publish_dir: 'test_paired_end'])

ch_multiqc_config = file("$launchDir/tests/data/multiqc/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

workflow test_multiqc {
    test_paired_end()

    input = [
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        test_paired_end.out.zip.collect { it[1] }.ifEmpty([])
    ]

    MULTIQC(*input)
}
