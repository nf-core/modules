#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALEVINFRY_GENERATEPERMITLIST } from '../../../../modules/alevinfry/generatepermitlist/main.nf'

workflow test_alevinfry_generatepermitlist {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    ALEVINFRY_GENERATEPERMITLIST ( input )
}
