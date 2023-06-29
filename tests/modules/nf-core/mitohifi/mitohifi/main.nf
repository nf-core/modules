#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI } from '../../../../../modules/nf-core/mitohifi/mitohifi/main.nf'

workflow test_mitohifi_mitohifi {
    
    species = "'Deilephila porcellus'"

    data_contigs = Channel.of([[id:"ilDeiPorc1"],[],file(params.test_data['eukaryotes']['deilephila_porcellus']['mito']['ilDeiPorc1.contigs.fa'], checkIfExists: true)])
    ref_gb = file(params.test_data['eukaryotes']['deilephila_porcells']['mito']['MW539688.1.gb'], checkIfExists: true)
    ref_gb = file(params.test_data['eukaryotes']['deilephila_porcells']['mito']['MW539688.1.fasta'], checkIfExists: true)
    code = 5
    MITOHIFI_MITOHIFI ( data_contigs, ref_fa, ref_gb, code )
}
