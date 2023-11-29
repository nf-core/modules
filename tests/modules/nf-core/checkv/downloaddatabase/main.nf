#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKV_DOWNLOADDATABASE } from '../../../../../modules/nf-core/checkv/downloaddatabase/main.nf'

workflow test_checkv_downloaddatabase {

    CHECKV_DOWNLOADDATABASE ()

}
