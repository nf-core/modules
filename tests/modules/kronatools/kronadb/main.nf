#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONATOOLS_KRONADB } from '../../../../modules/kronatools/kronadb/main.nf' addParams( options: [:] )

workflow test_kronatools_kronadb {
    KRONATOOLS_KRONADB ( )
}
