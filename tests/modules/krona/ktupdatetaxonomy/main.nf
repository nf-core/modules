#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONA_KTUPDATETAXONOMY } from '../../../../modules/krona/ktupdatetaxonomy/main.nf'

workflow test_krona_ktupdatetaxonomy {
    KRONA_KTUPDATETAXONOMY ( )
}
