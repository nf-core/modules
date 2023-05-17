#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../../../../modules/nf-core/deeptools/multibamsummary/main.nf'
include { DEEPTOOLS_PLOTCORRELATION } from '../../../../../modules/nf-core/deeptools/plotcorrelation/main.nf'

workflow test_plotcorrelation {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)],
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)],
        [ "test_bam1", "test_bam2" ]
    ]

    DEEPTOOLS_MULTIBAMSUMMARY ( input ) 

    correlation_method = 'spearman'
    correlation_plot_type = 'heatmap'
    DEEPTOOLS_PLOTCORRELATION ( DEEPTOOLS_MULTIBAMSUMMARY.out.matrix, correlation_method, correlation_plot_type )
}

workflow test_plotcorrelation_no_method {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)],
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)],
        [ "test_bam1", "test_bam2" ]
    ]

    DEEPTOOLS_MULTIBAMSUMMARY ( input ) 

    correlation_method = ''
    correlation_plot_type = 'heatmap'
    DEEPTOOLS_PLOTCORRELATION ( DEEPTOOLS_MULTIBAMSUMMARY.out.matrix, correlation_method, correlation_plot_type )
}

workflow test_plotcorrelation_no_plot_type {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)],
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)],
        [ "test_bam1", "test_bam2" ]
    ]

    DEEPTOOLS_MULTIBAMSUMMARY ( input ) 

    correlation_method = 'spearman'
    correlation_plot_type = ''
    DEEPTOOLS_PLOTCORRELATION ( DEEPTOOLS_MULTIBAMSUMMARY.out.matrix, correlation_method, correlation_plot_type )
}