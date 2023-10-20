#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPENMS_PEAKPICKERHIRES } from '../../../../../modules/nf-core/openms/peakpickerhires/main.nf'

workflow test_openms_peakpickerhires {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    OPENMS_PEAKPICKERHIRES ( input )
}
