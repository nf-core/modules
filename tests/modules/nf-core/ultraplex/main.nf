#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRAPLEX } from '../../../../modules/nf-core/ultraplex/main.nf'

workflow test_ultraplex {

    barcodes = file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/barcodes.csv", checkIfExists: true)
    input = [[id: "test"], file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/reads/multiplexed.fastq.gz", checkIfExists: true)]

    ULTRAPLEX (
        input,
        barcodes,
        ""
    )
}

workflow test_ultraplex_adaptor {

    barcodes = file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/barcodes.csv", checkIfExists: true)
    input = [[id: "test"], file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/reads/multiplexed.fastq.gz", checkIfExists: true)]

    ULTRAPLEX (
        input,
        barcodes,
        "AGATCGGAAGAGCGGTTCAG"
    )
}

workflow test_ultraplex_pairedend {

    barcodes = file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/barcodes.csv", checkIfExists: true)
    input = [[id: "test"], [file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/reads/multiplexed.fastq.gz", checkIfExists: true), file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/reads/multiplexed2.fastq.gz", checkIfExists: true)]]

    ULTRAPLEX (
        input,
        barcodes,
        ""
    )
}
