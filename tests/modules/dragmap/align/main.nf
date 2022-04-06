#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DRAGMAP_HASHTABLE } from '../../../../modules/dragmap/hashtable/main.nf'
include { DRAGMAP_ALIGN     } from '../../../../modules/dragmap/align/main.nf'

workflow test_dragmap_align_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    DRAGMAP_HASHTABLE ( fasta )
    DRAGMAP_ALIGN ( input, DRAGMAP_HASHTABLE.out.hashmap, false )
}

workflow test_dragmap_align_single_end_sort {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    DRAGMAP_HASHTABLE ( fasta )
    DRAGMAP_ALIGN ( input, DRAGMAP_HASHTABLE.out.hashmap, true )
}

workflow test_dragmap_align_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    DRAGMAP_HASHTABLE ( fasta )
    DRAGMAP_ALIGN ( input, DRAGMAP_HASHTABLE.out.hashmap, false )
}

workflow test_dragmap_align_paired_end_sort {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    DRAGMAP_HASHTABLE ( fasta )
    DRAGMAP_ALIGN ( input, DRAGMAP_HASHTABLE.out.hashmap, true )
}
