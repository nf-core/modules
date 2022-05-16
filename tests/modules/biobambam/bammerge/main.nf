#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BIOBAMBAM_BAMMERGE } from '../../../../modules/biobambam/bammerge/main.nf'

workflow test_biobambam_bammerge {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
        ]
    ]

    BIOBAMBAM_BAMMERGE ( input )
}
