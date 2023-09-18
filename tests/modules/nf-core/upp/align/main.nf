#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UPP_ALIGN } from '../../../../../modules/nf-core/upp/align/main.nf'

workflow test_upp_align {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    UPP_ALIGN ( input, [ [:], [ ] ] )
}
