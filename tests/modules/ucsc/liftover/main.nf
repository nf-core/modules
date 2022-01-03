#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_LIFTOVER } from '../../../../modules/ucsc/liftover/main.nf'

workflow test_ucsc_liftover {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)]
    chain =  file(params.test_data['homo_sapiens']['genome']['genome_chain_gz'], checkIfExists: true)

    UCSC_LIFTOVER ( input, chain )
}
