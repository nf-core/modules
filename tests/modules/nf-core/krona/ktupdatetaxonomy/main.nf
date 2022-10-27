#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { KRONA_KTUPDATETAXONOMY } from "$moduleDir/modules/nf-core/krona/ktupdatetaxonomy/main.nf"

workflow test_krona_ktupdatetaxonomy {
    KRONA_KTUPDATETAXONOMY ( )
}
