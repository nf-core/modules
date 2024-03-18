#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_COLLATEFASTQ } from '../../../../../modules/nf-core/samtools/collatefastq/main.nf'

workflow test_samtools_collatefastq_bam {
    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    fasta = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
            ]

    SAMTOOLS_COLLATEFASTQ ( input, fasta, false )
}

workflow test_samtools_collatefastq_bam_single_end {
    input = [
                [ id:'test', single_end:true ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    fasta = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
            ]

    SAMTOOLS_COLLATEFASTQ ( input, fasta, false )
}

workflow test_samtools_collatefastq_cram {
    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    fasta = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
            ]

    SAMTOOLS_COLLATEFASTQ ( input, fasta, false )
}

workflow test_samtools_collatefastq_bam_interleaved {
    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    fasta = [
                [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
            ]

    SAMTOOLS_COLLATEFASTQ ( input, fasta, true )
}
