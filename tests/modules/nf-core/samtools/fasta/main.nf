#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FASTA } from '../../../../../modules/nf-core/samtools/fasta/main.nf'

workflow test_samtools_fasta {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SAMTOOLS_FASTA ( input, false )
}

workflow test_samtools_fasta_interleave {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SAMTOOLS_FASTA ( input, true )
}
