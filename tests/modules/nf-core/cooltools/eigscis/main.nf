#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLTOOLS_EIGSCIS } from '../../../../../modules/nf-core/cooltools/eigscis/main.nf'

workflow test_cooltools_eigscis {

    input = [
        [ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/open2c/cooltools/master/tests/data/CN.mm9.1000kb.cool', checkIfExists: true),
        file('s3://ngi-igenomes/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa', checkIfExists: true),
        file('https://raw.githubusercontent.com/open2c/cooltools/master/tests/data/mm9.chrom.sizes.reduced', checkIfExists: true)
    ]

    COOLTOOLS_EIGSCIS ( input,  1000000)
}
