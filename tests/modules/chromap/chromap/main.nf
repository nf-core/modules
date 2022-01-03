#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHROMAP_INDEX                           } from '../../../../modules/chromap/index/main.nf'
include { CHROMAP_CHROMAP as CHROMAP_CHROMAP_BASE } from '../../../../modules/chromap/chromap/main.nf'
include { CHROMAP_CHROMAP as CHROMAP_CHROMAP_SAM  } from '../../../../modules/chromap/chromap/main.nf'

workflow test_chromap_chromap_single_end {

    // Test single-end and gz compressed output
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    CHROMAP_INDEX ( fasta )
    CHROMAP_CHROMAP_BASE (
        input,                      // meta + read data
        fasta,                      // reference genome
        CHROMAP_INDEX.out.index,    // reference index
        [],                         // barcode file
        [],                         // barcode whitelist
        [],                         // chromosome order file
        []                          // pairs chromosome order file
    )
}

workflow test_chromap_chromap_paired_end {

    // Test paired-end and gz compressed output
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    CHROMAP_INDEX ( fasta )
    CHROMAP_CHROMAP_BASE (
        input,                      // meta + read data
        fasta,                      // reference genome
        CHROMAP_INDEX.out.index,    // reference index
        [],                         // barcode file
        [],                         // barcode whitelist
        [],                         // chromosome order file
        []                          // pairs chromosome order file
    )
}

workflow test_chromap_chromap_paired_bam {

    // Test paired-end and bam output
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    CHROMAP_INDEX ( fasta )
    CHROMAP_CHROMAP_SAM (
        input,                      // meta + read data
        fasta,                      // reference genome
        CHROMAP_INDEX.out.index,    // reference index
        [],                         // barcode file
        [],                         // barcode whitelist
        [],                         // chromosome order file
        []                          // pairs chromosome order file
    )
}
