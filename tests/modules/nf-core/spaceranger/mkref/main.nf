#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPACERANGER_MKREF } from '../../../../../modules/nf-core/spaceranger/mkref/main.nf'

workflow test_spaceranger_mkref {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    SPACERANGER_MKREF ( input )
}
