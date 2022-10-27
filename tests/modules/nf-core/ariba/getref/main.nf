#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARIBA_GETREF } from "$moduleDir/modules/nf-core/ariba/getref/main.nf"

workflow test_ariba_getref {
    ARIBA_GETREF ( "card" )
}
