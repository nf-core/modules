#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_MATRIXFILTER } from '../../../../../modules/nf-core/custom/matrixfilter/main.nf'

empty_samplesheet = [[],[]]

workflow test_custom_matrixfilter {

    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)

    ch_samples_matrix = [ [ "id":"SRP254919" ], expression_matrix_file ]
    
    CUSTOM_MATRIXFILTER(
        ch_samples_matrix,
        empty_samplesheet
    )
}

workflow test_custom_matrixfilter_prop {

    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)

    ch_samples_matrix = [ [ "id":"SRP254919" ], expression_matrix_file ]
    
    CUSTOM_MATRIXFILTER(
        ch_samples_matrix,
        empty_samplesheet
    )
}

workflow test_custom_matrixfilter_group {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)

    ch_samples_matrix = [ [ "id":"SRP254919" ], expression_matrix_file ]
    ch_samplesheet = [ [ "id":"SRP254919" ], expression_sample_sheet ]    

    CUSTOM_MATRIXFILTER(
        ch_samples_matrix,
        ch_samplesheet
    )
}

workflow test_custom_matrixfilter_na_prop {

    expression_sample_sheet = file(params.test_data['proteomics']['maxquant']['mq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['proteomics']['maxquant']['mq_proteus_mat'], checkIfExists: true)

    ch_samples_matrix = [ [ "id":"mq_prop" ], expression_matrix_file ]
    ch_samplesheet = [ [ "id":"mq_prop" ], expression_sample_sheet ]

    CUSTOM_MATRIXFILTER(
        ch_samples_matrix,
        ch_samplesheet
    )
}

workflow test_custom_matrixfilter_na_samples {

    expression_sample_sheet = file(params.test_data['proteomics']['maxquant']['mq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['proteomics']['maxquant']['mq_proteus_mat'], checkIfExists: true)

    ch_samples_matrix = [ [ "id":"mq_samples" ], expression_matrix_file ]
    ch_samplesheet = [ [ "id":"mq_samples" ], expression_sample_sheet ]

    CUSTOM_MATRIXFILTER(
        ch_samples_matrix,
        ch_samplesheet
    )
}

workflow test_custom_matrixfilter_var {

    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)

    ch_samples_matrix = [ [ "id":"SRP254919" ], expression_matrix_file ]
    
    CUSTOM_MATRIXFILTER(
        ch_samples_matrix,
        empty_samplesheet
    )
}
