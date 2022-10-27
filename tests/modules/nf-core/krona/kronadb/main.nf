#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { KRONA_KRONADB } from "$moduleDir/modules/nf-core/krona/kronadb/main.nf"

workflow test_krona_kronadb {
    KRONA_KRONADB ( )
}
