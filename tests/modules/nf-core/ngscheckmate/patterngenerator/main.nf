#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_BUILD } from '../../../../../modules/nf-core/bowtie/build/main.nf'
include { NGSCHECKMATE_PATTERNGENERATOR } from '../../../../../modules/nf-core/ngscheckmate/patterngenerator/main.nf'

workflow test_ngscheckmate_patterngenerator {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    fasta    = [
        [ id: 'sarscov2' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]

    bowtie_fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BOWTIE_BUILD ( bowtie_fasta )

    NGSCHECKMATE_PATTERNGENERATOR ( input, fasta, BOWTIE_BUILD.out.index )
}
