#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LIMMA_DIFFERENTIAL } from '../../../../../modules/nf-core/limma/differential/main.nf'

workflow test_limma_differential {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)

    Channel.fromPath(expression_contrasts)
        .splitCsv ( header:true, sep:',' )
        .map{
            tuple(it, expression_sample_sheet, expression_matrix_file)
        }
        .set{
            input
        }

    LIMMA_DIFFERENTIAL (
        input
    )
}

