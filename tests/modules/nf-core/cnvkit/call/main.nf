#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_CALL } from '../../../../../modules/nf-core/cnvkit/call/main.nf'

workflow test_cnvkit_call {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cns'], checkIfExists: true),
       []
    ]

    CNVKIT_CALL ( input )
}

workflow test_cnvkit_call_with_vcf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cns'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true)
    ]

    CNVKIT_CALL ( input )
}