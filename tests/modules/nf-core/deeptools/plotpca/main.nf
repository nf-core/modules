#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../../../../modules/nf-core/deeptools/multibamsummary/main.nf'
include { DEEPTOOLS_PLOTPCA         } from '../../../../../modules/nf-core/deeptools/plotpca/main.nf'

workflow test_plotpca {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)],
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)],
        [ "test_bam1", "test_bam2" ]
    ]

    DEEPTOOLS_MULTIBAMSUMMARY ( input )
    DEEPTOOLS_PLOTPCA ( DEEPTOOLS_MULTIBAMSUMMARY.out.matrix )
}
