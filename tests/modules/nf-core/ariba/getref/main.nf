#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARIBA_GETREF } from '../../../../modules/ariba/getref/main.nf'

workflow test_ariba_getref {
    ARIBA_GETREF ( "card" )
}
