#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_JACCARD } from '../../../../../modules/nf-core/bedtools/jaccard/main.nf'

workflow test_bedtools_jaccard {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true)
    ]

    BEDTOOLS_JACCARD ( input, [[],[]] )
}

workflow test_bedtools_jaccard_genome {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true)
    ]

    genome = [
        [ id:'genome' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    BEDTOOLS_JACCARD ( input, genome )
}
