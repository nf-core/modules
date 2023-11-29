#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_TARGET } from '../../../../../modules/nf-core/cnvkit/target/main.nf'

workflow test_cnvkit_target {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]
    annotations = [
        [:],
        []
    ]

    CNVKIT_TARGET ( input, annotations )
}

workflow test_cnvkit_target_with_gff3 {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]
    annotations = [
        [ id: "annotations" ],
        file(params.test_data['homo_sapiens']['genome']['genome_gff3'], checkIfExists: true)
    ]

    CNVKIT_TARGET ( input, annotations )
}
