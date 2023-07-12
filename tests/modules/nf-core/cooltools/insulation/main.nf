#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLTOOLS_INSULATION } from '../../../../../modules/nf-core/cooltools/insulation/main.nf'

workflow test_cooltools_insulation {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    COOLTOOLS_INSULATION ( input )
}
