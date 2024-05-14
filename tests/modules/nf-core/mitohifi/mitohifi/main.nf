#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_CONTIGS } from '../../../../../modules/nf-core/mitohifi/mitohifi/main.nf'
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_READS } from '../../../../../modules/nf-core/mitohifi/mitohifi/main.nf'

workflow test_mitohifi_mitohifi_contigs {
    
    species = "'Deilephila porcellus'"

    data_contigs = Channel.of([[id:"ilDeiPorc1"],file(params.test_data['deilephila_porcellus']['mito']['contigs'], checkIfExists: true)])
    ref_gb = file(params.test_data['deilephila_porcellus']['mito']['ref_gb'], checkIfExists: true)
    ref_fa = file(params.test_data['deilephila_porcellus']['mito']['ref_fa'], checkIfExists: true)
    input_mode = Channel.of('c')
    code = 5
    MITOHIFI_MITOHIFI_CONTIGS ( data_contigs, ref_fa, ref_gb, input_mode, code )
}


workflow test_mitohifi_mitohifi_reads {

    species = "'Deilephila porcellus'"

    data_reads = Channel.of([[id:"ilDeiPorc1"],file(params.test_data['deilephila_porcellus']['mito']['hifi_reads'], checkIfExists: true)])
    ref_gb = file(params.test_data['deilephila_porcellus']['mito']['ref_gb'], checkIfExists: true)
    ref_fa = file(params.test_data['deilephila_porcellus']['mito']['ref_fa'], checkIfExists: true)
    input_mode = Channel.of('r')
    code = 5
    MITOHIFI_MITOHIFI_READS ( data_reads, ref_fa, ref_gb, input_mode, code )
}
