#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_SCATTERINTERVALSBYNS } from '../../../../../modules/nf-core/picard/scatterintervalsbyns/main.nf'

workflow test_picard_scatterintervalsbyns {
    
    fasta = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fai = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    dict = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    ]

    PICARD_SCATTERINTERVALSBYNS ( fasta, fai, dict )
}
