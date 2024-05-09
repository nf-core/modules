#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PILON } from '../../../../modules/nf-core/pilon/main.nf'

workflow test_pilon {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    bam_tuple_ch = Channel.of([ [ id:'test', single_end:false ], // meta map
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                                ])

    PILON ( input, bam_tuple_ch, "bam" )
}
