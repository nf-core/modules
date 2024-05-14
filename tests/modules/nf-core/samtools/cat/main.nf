#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_CAT } from '../../../../../modules/nf-core/samtools/cat/main.nf'

workflow test_samtools_cat {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_unaligned_bam'], checkIfExists: true)
    ]

    SAMTOOLS_CAT ( input )
}
