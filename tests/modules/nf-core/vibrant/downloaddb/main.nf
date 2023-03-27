#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VIBRANT_DOWNLOADDB } from '../../../../../modules/nf-core/vibrant/downloaddb/main.nf'

workflow test_vibrant_downloaddb {
    VIBRANT_DOWNLOADDB ( )
}
