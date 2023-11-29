#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_COVERAGE } from '../../../../../modules/nf-core/bedtools/coverage/main.nf'

workflow test_bedtools_coverage {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    BEDTOOLS_COVERAGE ( input, [] )
}

workflow test_bedtools_coverage_fasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    fasta_fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    BEDTOOLS_COVERAGE ( input, fasta_fai )
}
