#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPEFF } from '../../../modules/snpeff/main.nf'

workflow test_snpeff {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    SNPEFF ( input, "WBcel235.99", [] )
}
