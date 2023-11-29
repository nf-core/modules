#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAKTA_BAKTADBDOWNLOAD } from '../../../../../modules/nf-core/bakta/baktadbdownload/main.nf'

workflow test_bakta_baktadbdownload {
    BAKTA_BAKTADBDOWNLOAD ( )
}
