#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SHIFTFASTA               } from '../../../../../modules/nf-core/gatk4/shiftfasta/main.nf'

workflow test_gatk4_shiftfasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    index =  file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict  =  file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_SHIFTFASTA ( input, index, dict )
}
