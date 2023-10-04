#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROPR_PROPD } from '../../../../../modules/nf-core/propr/propd/main.nf'

workflow test_propr_propd_default {
    
    matrix = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    ]
    samplesheet = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    ]

    PROPR_PROPD ( matrix, samplesheet )
}

workflow test_propr_propd_default_permutation {
    
    matrix = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    ]
    samplesheet = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    ]

    PROPR_PROPD ( matrix, samplesheet )
}

workflow test_propr_propd_default_boxcox_permutation {
    
    matrix = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    ]
    samplesheet = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    ]

    PROPR_PROPD ( matrix, samplesheet )
}

workflow test_propr_propd_thetae_permutation {
    
    matrix = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    ]
    samplesheet = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    ]

    PROPR_PROPD ( matrix, samplesheet )
}

workflow test_propr_propd_thetae_boxcox_permutation {
    
    matrix = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    ]
    samplesheet = [
        [ id:'test' ],
        file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    ]

    PROPR_PROPD ( matrix, samplesheet )
}