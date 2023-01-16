#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_ANCESTRY } from '../../../../../modules/nf-core/somalier/ancestry/main.nf'

workflow test_somalier_ancestry {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SOMALIER_ANCESTRY ( input )
}
