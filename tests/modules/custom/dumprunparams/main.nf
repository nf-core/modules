#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_DUMPRUNPARAMS } from '../../../../modules/custom/dumprunparams/main.nf'
include { MULTIQC } from '../../../modules/multiqc/main.nf'


workflow test_custom_dumprunparams_multi {
    CUSTOM_DUMPRUNPARAMS ( [] )
    MULTIQC ( CUSTOM_DUMPRUNPARAMS.out )
}


workflow test_custom_dumprunparams_oneexclude {
    CUSTOM_DUMPRUNPARAMS ( [ 'test_data' ] )
    MULTIQC ( CUSTOM_DUMPRUNPARAMS.out )
}

workflow test_custom_dumprunparams_twoexclude {
    CUSTOM_DUMPRUNPARAMS ( [ 'test_data', 'enable_conda' ] )
    MULTIQC ( CUSTOM_DUMPRUNPARAMS.out )
}
