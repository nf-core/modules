#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKV_UPDATEDATABASE } from '../../../../../modules/nf-core/checkv/updatedatabase/main.nf'

workflow test_checkv_updatedatabase {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    CHECKV_UPDATEDATABASE ( input )
}
