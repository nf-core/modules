#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_MPILEUP } from '../../../../modules/samtools/mpileup/main.nf'

workflow test_samtools_mpileup {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                []
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    SAMTOOLS_MPILEUP ( input, fasta )
}

workflow test_samtools_mpileup_intervals {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    SAMTOOLS_MPILEUP ( input, fasta )
}
