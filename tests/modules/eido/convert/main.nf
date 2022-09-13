#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EIDO_CONVERT } from '../../../../modules/eido/convert/main.nf'

workflow test_eido_convert_nextflow_samplesheet {

    nextflow_samplesheet = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/pep/test_nextflow_original_samplesheet.csv", checkIfExists: true)

    EIDO_CONVERT ( nextflow_samplesheet )
}


workflow test_eido_convert_pep_project {

    pep_project = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/pep/test_pep_format_files/config.yaml", checkIfExists: true)

    EIDO_CONVERT ( pep_project )
}
