#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_COLLATE } from '../../../../../modules/nf-core/samtools/collate/main.nf'

workflow test_samtools_collate {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SAMTOOLS_COLLATE ( input, [] )
}

workflow test_samtools_collate_cram {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    SAMTOOLS_COLLATE ( input, fasta )
}
