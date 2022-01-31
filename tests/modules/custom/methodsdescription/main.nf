#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_METHODSDESCRIPTION } from '../../../../modules/custom/methodsdescription/main.nf'
include { MULTIQC                   } from '../../../../modules/multiqc/main.nf'


workflow test_custom_methodsdescription {

    CUSTOM_METHODSDESCRIPTION ( )
    MULTIQC( CUSTOM_METHODSDESCRIPTION.out.mqc_html )

}
