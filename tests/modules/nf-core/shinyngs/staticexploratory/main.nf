#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHINYNGS_STATICEXPLORATORY } from '../../../../../modules/nf-core/shinyngs/staticexploratory/main.nf'

workflow test_shinyngs_staticexploratory {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true) 
    expression_feature_meta = file(params.test_data['mus_musculus']['genome']['rnaseq_genemeta'], checkIfExists: true) 
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true) 

    SHINYNGS_STATICEXPLORATORY ( 
        [ [ "id":"treatment" ], expression_sample_sheet, expression_feature_meta, [ expression_matrix_file ] ],
    )
}

workflow test_shinyngs_staticexploratory_html {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true) 
    expression_feature_meta = file(params.test_data['mus_musculus']['genome']['rnaseq_genemeta'], checkIfExists: true) 
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true) 

    SHINYNGS_STATICEXPLORATORY ( 
        [ [ "id":"treatment" ], expression_sample_sheet, expression_feature_meta, [ expression_matrix_file ] ],
    )
}
