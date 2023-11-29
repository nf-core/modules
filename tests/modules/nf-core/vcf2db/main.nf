#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF2DB } from '../../../../modules/nf-core/vcf2db/main.nf'

workflow test_vcf2db {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['justhusky_minimal_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['justhusky_ped'], checkIfExists: true)
    ]

    VCF2DB ( input )
}
