#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_DEPTH } from '../../../modules/samtools/depth/main.nf'
include { SEXDETERRMINE } from '../../../modules/sexdeterrmine/main.nf'

workflow test_sexdeterrmine {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/eager/testdata/Human/bam/JK2067_downsampled_s0.1.bam", checkIfExists: true) ]

    SAMTOOLS_DEPTH ( input )
    SEXDETERRMINE ( SAMTOOLS_DEPTH.out.tsv, [] )
}
