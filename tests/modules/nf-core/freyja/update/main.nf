#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREYJA_UPDATE } from '../../../../../modules/nf-core/freyja/update/main.nf'

workflow test_freyja_update {
    FREYJA_UPDATE ()
}
