#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INITIALISE } from '../../../../subworkflows/nf-core/initialise/main.nf'

workflow test_initialise {
    INITIALISE ( 
        params.version,
        params.help,
        params.validate_params
    )
}
