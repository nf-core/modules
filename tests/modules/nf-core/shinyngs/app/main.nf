#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHINYNGS_APP } from '../../../../../modules/nf-core/shinyngs/app/main.nf'

workflow test_shinyngs_app {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true) 
    expression_feature_meta = file(params.test_data['mus_musculus']['genome']['rnaseq_genemeta'], checkIfExists: true) 
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true) 
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true) 
    expression_differential = file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true) 

    SHINYNGS_APP ( 
        expression_sample_sheet,
        expression_feature_meta,
        [ expression_matrix_file ],
        expression_contrasts,
        [ expression_differential ] 
    )
}
