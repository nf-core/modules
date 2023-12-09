#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IPHOP_DOWNLOAD } from '../../../../../modules/nf-core/iphop/download/main.nf'

workflow test_iphop_download {

    IPHOP_DOWNLOAD ( )

}
