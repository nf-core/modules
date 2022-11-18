#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MATRIXFILTER } from '../../../../modules/nf-core/matrixfilter/main.nf'

workflow test_matrixfilter {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)

    ch_samples_matrix = [ [ "id":"SRP254919" ], expression_sample_sheet, expression_matrix_file ]
    
    MATRIXFILTER(
        ch_samples_matrix
    )
}

workflow test_matrixfilter_prop {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)

    ch_samples_matrix = [ [ "id":"SRP254919" ], expression_sample_sheet, expression_matrix_file ]
    
    MATRIXFILTER(
        ch_samples_matrix
    )
}

workflow test_matrixfilter_group {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)

    ch_samples_matrix = [ [ "id":"SRP254919" ], expression_sample_sheet, expression_matrix_file ]
    
    MATRIXFILTER(
        ch_samples_matrix
    )
}
