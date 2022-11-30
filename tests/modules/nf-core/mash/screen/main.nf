#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MASH_SKETCH } from '../../../../../modules/nf-core/mash/sketch/main.nf'
include { MASH_SCREEN } from '../../../../../modules/nf-core/mash/screen/main.nf'

workflow test_mash_screen {

    input = [
        [ id:'test', single_end:false], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    sars_db = [
        [ id: 'sars_db' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    MASH_SKETCH ( sars_db )
    MASH_SCREEN ( input, MASH_SKETCH.out.mash.map { meta, sketch -> sketch } )
}
