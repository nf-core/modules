#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MSISENSORPRO_SCAN } from '../../../../modules/msisensorpro/scan/main.nf'

workflow test_msisensorpro_scan {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    MSISENSORPRO_SCAN ( input )
}
