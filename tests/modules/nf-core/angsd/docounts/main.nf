#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANGSD_DOCOUNTS } from '../../../../../modules/nf-core/angsd/docounts/main.nf'

workflow test_angsd_docounts {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam_bai'], checkIfExists: true),
        []
    ]

    ANGSD_DOCOUNTS ( input )
}
