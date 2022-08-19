#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VERIFYBAMID } from '../../../modules/verifybamid/main.nf'

workflow test_verifybamid {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    VERIFYBAMID ( input )
}
