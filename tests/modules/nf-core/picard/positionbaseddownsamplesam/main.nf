#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_POSITIONBASEDDOWNSAMPLESAM } from '../../../../../modules/nf-core/picard/positionbaseddownsamplesam/main.nf'

workflow test_picard_positionbaseddownsamplesam {

    input = [
        [ id:'test', single_end:false], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        100,
        0.5
    ]

    PICARD_POSITIONBASEDDOWNSAMPLESAM ( input )
}
