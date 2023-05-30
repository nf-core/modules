#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_PREPROCESS } from '../../../../../modules/nf-core/varlociraptor/preprocess/main.nf'

workflow test_varlociraptor_preprocess {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        file("./test.alignment-properties.json")
    ]

    fasta = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fai= [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_PREPROCESS ( input, fasta, fai )
}
