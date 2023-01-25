#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHINYNGS_VALIDATEFOMCOMPONENTS } from '../../../../../modules/nf-core/shinyngs/validatefomcomponents/main.nf'

workflow test_shinyngs_validatefomcomponents {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_feature_meta = file(params.test_data['mus_musculus']['genome']['rnaseq_genemeta'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)

    SHINYNGS_VALIDATEFOMCOMPONENTS (
        [ [ "id":"SRP254919" ], expression_sample_sheet, expression_matrix_file ],
        [ [ "id":"SRP254919" ], expression_feature_meta ],
        [ [ "id":"SRP254919" ], expression_contrasts ]
    )
}

workflow test_shinyngs_validatefomcomponents_no_feature {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)

    SHINYNGS_VALIDATEFOMCOMPONENTS (
        [ [ "id":"SRP254919" ], expression_sample_sheet, expression_matrix_file ],
        [ [ "id":"SRP254919" ], [] ],
        [ [ "id":"SRP254919" ], expression_contrasts ]
    )
}
