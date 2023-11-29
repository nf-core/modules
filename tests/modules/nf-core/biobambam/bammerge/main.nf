#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BIOBAMBAM_BAMMERGE } from '../../../../../modules/nf-core/biobambam/bammerge/main.nf'

workflow test_biobambam_bammerge_paired {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
        ]
    ]

    BIOBAMBAM_BAMMERGE ( input )
}

workflow test_biobambam_bammerge_single {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
        ]
    ]

    BIOBAMBAM_BAMMERGE ( input )
}
