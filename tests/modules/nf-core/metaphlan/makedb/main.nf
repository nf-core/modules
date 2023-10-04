#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METAPHLAN_MAKEDB } from '../../../../../modules/nf-core/metaphlan/makedb/main.nf'

workflow test_metaphlan_makedb {

    METAPHLAN_MAKEDB ()
}
