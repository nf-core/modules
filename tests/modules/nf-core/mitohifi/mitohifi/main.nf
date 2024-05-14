#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI } from '../../../../../modules/nf-core/mitohifi/mitohifi/main.nf'

workflow test_mitohifi_mitohifi {
    
    species = "'Deilephila porcellus'"

    data_contigs = Channel.of([[id:"ilDeiPorc1"],[],file(params.test_data['deilephila_porcellus']['mito']['contigs'], checkIfExists: true)])
    ref_gb = file(params.test_data['deilephila_porcellus']['mito']['ref_gb'], checkIfExists: true)
    ref_fa = file(params.test_data['deilephila_porcellus']['mito']['ref_fa'], checkIfExists: true)
    code = 5
    MITOHIFI_MITOHIFI ( data_contigs, ref_fa, ref_gb, code )
}
