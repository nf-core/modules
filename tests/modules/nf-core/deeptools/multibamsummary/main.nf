#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../../../../modules/nf-core/deeptools/multibamsummary/main.nf'

workflow test_multibamsummary_bam {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        [ "test_bam" ]
    ]

    DEEPTOOLS_MULTIBAMSUMMARY ( input ) 
}

workflow test_multibamsummary_bam_no_label {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        [ ] // error
    ]

    DEEPTOOLS_MULTIBAMSUMMARY ( input ) 
}