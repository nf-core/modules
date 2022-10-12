#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SMOOVE } from '../../../../modules/nf-core/smoove/main.nf'

workflow test_smoove {
    
    bam_tuple_ch = Channel.of([ [ id:'test' ], // meta map
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
<<<<<<< HEAD
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
=======
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
>>>>>>> 56aa4e3c475b756088625940eb5f5f631544f2e7
    				])
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    SMOOVE ( bam_tuple_ch, fasta, fai )
}
