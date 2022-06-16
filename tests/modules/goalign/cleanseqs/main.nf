#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GOALIGN_CLEANSEQS } from '../../../../modules/goalign/cleanseqs/main.nf'

workflow test_goalign_cleanseqs {
    
    input = file(params.test_data['sarscov2']['genome']['all_sites_fas'], checkIfExists: true)

    GOALIGN_CLEANSEQS ( input )
}
