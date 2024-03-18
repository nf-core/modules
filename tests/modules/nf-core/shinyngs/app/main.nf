#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHINYNGS_APP } from '../../../../../modules/nf-core/shinyngs/app/main.nf'

workflow test_shinyngs_app_multi_matrix {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true) 
    expression_feature_meta = file(params.test_data['mus_musculus']['genome']['rnaseq_genemeta'], checkIfExists: true) 
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true) 

    // Load the same matrix twice to immitate multiple processings (e.g. raw
    // and normalised)
    
    raw_expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true) 
    raw_expression_matrix_file.copyTo('normalised.tsv') 
    normalised_expression_matrix_file = file('normalised.tsv')
    
    // Same for differential stats - we need two tables for the (now) two contrasts 

    expression_differential = file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true) 
    expression_differential.copyTo('second_contrast_stats.tsv') 
    second_contrast_stats = file('second_contrast_stats.tsv')
    contrast_stats_assay = Channel.value(1)

    SHINYNGS_APP ( 
        [ [ "id":"SRP254919" ], expression_sample_sheet,  expression_feature_meta, [ raw_expression_matrix_file, normalised_expression_matrix_file ] ],
        [ [ "id":"SRP254919" ], expression_contrasts, [ expression_differential, second_contrast_stats  ] ],
        contrast_stats_assay
    )
}

workflow test_shinyngs_app_single_matrix {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true) 
    expression_feature_meta = file(params.test_data['mus_musculus']['genome']['rnaseq_genemeta'], checkIfExists: true) 
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true) 
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true) 
    expression_differential = file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true) 

    // We need two tables for the (now) two contrasts 

    expression_differential = file(params.test_data['mus_musculus']['genome']['deseq_results'], checkIfExists: true) 
    expression_differential.copyTo('second_contrast_stats.tsv') 
    second_contrast_stats = file('second_contrast_stats.tsv')
    contrast_stats_assay = Channel.value(1)

    SHINYNGS_APP ( 
        [ [ "id":"SRP254919" ], expression_sample_sheet, expression_feature_meta,  [ expression_matrix_file ] ],
        [ [ "id":"SRP254919" ], expression_contrasts, [ expression_differential, second_contrast_stats  ] ],
        contrast_stats_assay
    )
}

