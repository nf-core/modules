#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISCUIT_VCF2BED } from '../../../../modules/biscuit/vcf2bed/main.nf'

workflow test_biscuit_vcf2bed {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/biscuit/test.vcf.gz', checkIfExists: true)
    ]

    BISCUIT_VCF2BED ( input )

}
