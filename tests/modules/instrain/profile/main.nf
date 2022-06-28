#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INSTRAIN_PROFILE } from '../../../modules/instrain/profile/main.nf'

workflow test_instrain_profile {
    
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['bacteroides_fragilis']['illumina']['test1_paired_end_sorted_bam'], checkIfExists: true)
        ]
    ]
    genome_fasta = file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)

    INSTRAIN_PROFILE ( input )
}
