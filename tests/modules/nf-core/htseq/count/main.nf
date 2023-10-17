#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HTSEQ_COUNT } from '../../../../../modules/nf-core/htseq/count/main.nf'

workflow test_htseq_count {

    bam = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)

    ]
    gtf = [
        [ id:'test2'], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    ]


    HTSEQ_COUNT ( bam,gtf)
}
