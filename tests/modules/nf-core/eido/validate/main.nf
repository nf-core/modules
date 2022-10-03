#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EIDO_VALIDATE } from '../../../../../modules/nf-core/eido/validate/main.nf'

workflow test_eido_validate_on_nextflow_samplesheet {

    samplesheet = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/pep/test_nextflow_original_samplesheet.csv", checkIfExists: true)
    schema = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/pep/test_samplesheet_schema.yaml", checkIfExists: true)
    pep_input_base_dir = []

    EIDO_VALIDATE ( samplesheet, schema, pep_input_base_dir )
}

workflow test_eido_validate_on_pep_config {

    samplesheet = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/pep/test_pep_format_files/config.yaml", checkIfExists: true)
    schema = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/pep/test_samplesheet_schema.yaml", checkIfExists: true)
    pep_input_base_dir = []

    EIDO_VALIDATE ( samplesheet, schema, pep_input_base_dir )
}
