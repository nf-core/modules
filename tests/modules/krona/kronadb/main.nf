#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONA_KRONADB } from '../../../../modules/krona/kronadb/main.nf'

workflow test_krona_kronadb {
    KRONA_KRONADB ( )
}
