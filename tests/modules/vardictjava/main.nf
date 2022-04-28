#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARDICTJAVA } from '../../../modules/vardictjava/main.nf'

workflow test_vardictjava {
    
    bam_input_ch = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ])

    bed = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    VARDICTJAVA ( bam_input_ch, bed, fasta )
}
