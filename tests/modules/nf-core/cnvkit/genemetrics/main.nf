#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_GENEMETRICS } from '../../../../../modules/nf-core/cnvkit/genemetrics/main.nf'

workflow test_cnvkit_genemetrics_with_cns {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cnr'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cns'], checkIfExists: true)
    ]

    CNVKIT_GENEMETRICS ( input )
}

workflow test_cnvkit_genemetrics_without_cns {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['cnvkit']['amplicon_cnr'], checkIfExists: true),
        []
    ]

    CNVKIT_GENEMETRICS ( input )
}