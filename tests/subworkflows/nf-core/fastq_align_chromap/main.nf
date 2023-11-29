#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHROMAP_INDEX       } from '../../../../modules/nf-core/chromap/index/main.nf'
include { FASTQ_ALIGN_CHROMAP } from '../../../../subworkflows/nf-core/fastq_align_chromap/main.nf'

workflow test_fastq_align_chromap_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = Channel.of([
        [ id:'test' ],
        [
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]
    ])

    CHROMAP_INDEX ( fasta )
    FASTQ_ALIGN_CHROMAP ( input, CHROMAP_INDEX.out.index, fasta, [], [], [], [] )
}

workflow test_fastq_align_chromap_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = Channel.of([
        [ id:'test' ],
        [
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]
    ])

    CHROMAP_INDEX ( fasta )
    FASTQ_ALIGN_CHROMAP ( input, CHROMAP_INDEX.out.index, fasta, [], [], [], [] )
}
