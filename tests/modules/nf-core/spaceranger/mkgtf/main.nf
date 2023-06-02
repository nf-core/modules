#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPACERANGER_MKGTF } from '../../../../../modules/nf-core/spaceranger/mkgtf/main.nf'

workflow test_spaceranger_mkgtf {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    SPACERANGER_MKGTF ( input )
}
