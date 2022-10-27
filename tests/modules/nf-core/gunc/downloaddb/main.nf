#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNC_DOWNLOADDB } from "$moduleDir/modules/nf-core/gunc/downloaddb/main.nf"

workflow test_gunc_downloaddb {

    input = 'progenomes'

    GUNC_DOWNLOADDB ( input )
}
