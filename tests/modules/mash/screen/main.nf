#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MASH_SKETCH } from '../../../../modules/mash/sketch/main.nf'
include { MASH_SCREEN } from '../../../../modules/mash/screen/main.nf'

workflow test_mash_screen {

    input = [
        [ id:'test', single_end:false], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fastx_db = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    MASH_SKETCH ( input )
    MASH_SCREEN ( MASH_SKETCH.out.mash, fastx_db )
}
