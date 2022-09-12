#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EIDO_VALIDATE } from '../../../../modules/eido/validate/main.nf'

workflow test_eido_validate_on_nextflow_samplesheet {

    samplesheet = file(params.test_data['pep']['nextflow_samplesheet'], checkIfExists: true)
    schema = file(params.test_data['pep']['schema'], checkIfExists: true)

    EIDO_VALIDATE ( samplesheet, schema )
}

workflow test_eido_validate_on_pep_config {

    samplesheet = file(params.test_data['pep']['pep_config'], checkIfExists: true)
    schema = file(params.test_data['pep']['schema'], checkIfExists: true)

    EIDO_VALIDATE ( samplesheet, schema )
}
