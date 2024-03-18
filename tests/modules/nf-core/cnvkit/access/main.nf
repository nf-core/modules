#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_ACCESS } from '../../../../../modules/nf-core/cnvkit/access/main.nf'

workflow test_cnvkit_access {


    fasta = [
        [ id: "test" ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    exclude = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
    ]

    CNVKIT_ACCESS ( fasta, exclude )
}

workflow test_cnvkit_access_no_exclude {


    fasta = [
        [ id: "test" ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    exclude = [
        [ id:'test' ], // meta map
        []
    ]

    CNVKIT_ACCESS ( fasta, exclude )
}

workflow test_cnvkit_access_multiple_exclude {


    fasta = [
        [ id: "test" ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    exclude = [
        [ id:'test' ], // meta map
        [
            file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_antitarget_bed'], checkIfExists: true)
        ]
    ]

    CNVKIT_ACCESS ( fasta, exclude )
}
