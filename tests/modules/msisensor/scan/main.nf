#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MSISENSOR_SCAN } from '../../../../modules/msisensor/scan/main.nf'

workflow test_msisensor_scan {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]

    MSISENSOR_SCAN ( input )
}
