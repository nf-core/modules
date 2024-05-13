#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES } from '../../../../../modules/nf-core/varlociraptor/estimatealignmentproperties/main.nf'

workflow test_varlociraptor_estimatealignmentproperties {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    fasta = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]


    fai= [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES ( input, fasta, fai )
}
