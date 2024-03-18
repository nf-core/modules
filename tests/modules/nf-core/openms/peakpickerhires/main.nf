#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPENMS_PEAKPICKERHIRES } from '../../../../../modules/nf-core/openms/peakpickerhires/main.nf'

workflow test_openms_peakpickerhires {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['proteomics']['msspectra']['profile_spectra'], checkIfExists: true)
    ]

    OPENMS_PEAKPICKERHIRES ( input )
}
