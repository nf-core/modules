#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EIDO_CONVERT } from '../../../../../modules/nf-core/eido/convert/main.nf'

workflow test_eido_convert_nextflow_samplesheet {

    nextflow_samplesheet = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/pep/test_nextflow_original_samplesheet.csv", checkIfExists: true)
    format = "csv"
    pep_input_base_dir = []

    EIDO_CONVERT ( nextflow_samplesheet, format, pep_input_base_dir )
}


workflow test_eido_convert_pep_project {

    pep_project = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/pep/test_pep_format_files/config.yaml", checkIfExists: true)
    format = "csv"
    pep_input_base_dir = []

    EIDO_CONVERT ( pep_project, format, pep_input_base_dir )
}
