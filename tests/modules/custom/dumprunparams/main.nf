#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_DUMPRUNPARAMS } from '../../../../modules/custom/dumprunparams/main.nf'
include { MULTIQC } from '../../../modules/multiqc/main.nf'


workflow test_custom_dumprunparams {

    CUSTOM_DUMPRUNPARAMS ( )
    MULTIQC ( CUSTOM_DUMPRUNPARAMS.out )
}
