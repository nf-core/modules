#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UPD } from '../../../../modules/nf-core/upd/main.nf'

workflow test_upd {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['justhusky_minimal_vcf_gz'], checkIfExists: true)
    ]

    UPD ( input )
}
