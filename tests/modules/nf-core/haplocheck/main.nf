#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPLOCHECK } from '../../../../modules/nf-core/haplocheck/main.nf'

workflow test_haplocheck {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_mito_vcf'], checkIfExists: true)
    ]

    HAPLOCHECK ( input )
}
