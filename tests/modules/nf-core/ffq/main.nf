#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FFQ } from '../../../../modules/nf-core/ffq/main.nf'

workflow test_ffq_single_id {

    FFQ ( [ 'SRR9990627' ] )
}

workflow test_ffq_multiple_ids {

    FFQ ( [ 'SRR9990627', 'SRX7347523' ] )
}
