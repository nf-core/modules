#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRAHUMANSCRUBBER_INITDB } from '../../../../../modules/nf-core/srahumanscrubber/initdb/main.nf'

workflow test_srahumanscrubber_initdb {
    SRAHUMANSCRUBBER_INITDB ( )
}
