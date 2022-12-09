#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPTOOLS_BAMCOVERAGE } from '../../../../../modules/nf-core/deeptools/bamcoverage/main.nf'

workflow test_deeptools_bamcoverage_bam {

    input = [
        [ id:'test', single_end:false ], // meta map
          file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
        ]

    DEEPTOOLS_BAMCOVERAGE ( input, [], [] )
}

workflow test_deeptools_bamcoverage_cram {

    input = [
        [ id:'test', single_end:false ], // meta map
          file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
        ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    DEEPTOOLS_BAMCOVERAGE ( input, fasta, fasta_fai)
}

workflow test_deeptools_bamcoverage_cram_no_fasta_fai {

    input = [
        [ id:'test', single_end:false ], // meta map
          file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
        ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    DEEPTOOLS_BAMCOVERAGE ( input, fasta, [])
}
