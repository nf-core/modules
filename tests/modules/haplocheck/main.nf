#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPLOCHECK } from '../../../modules/haplocheck/main.nf'

workflow test_haplocheck {

    input = [
        [ id:'test' ], // meta map
        // TODO: change to "params.test_data[]"
        file("https://github.com/nf-core/test-datasets/raw/5101234ce48c3eb08adeed922e30a6e57e4fe5fb/testdata/NA12878_chrM.vcf.gz", checkIfExists: true)
    ]

    HAPLOCHECK ( input )
}
