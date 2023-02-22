#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPEFF_SNPEFF } from '../../../../../modules/nf-core/snpeff/snpeff/main.nf'

workflow test_snpeff_snpeff {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    SNPEFF_SNPEFF ( input, "WBcel235.105", [] )
}
