#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INSTRAIN_PROFILE } from '../../../../../modules/nf-core/instrain/profile/main.nf'

workflow test_instrain_profile {
    
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
        ]
    ]
    genome_fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    INSTRAIN_PROFILE ( input , genome_fasta , [] , [] )
}
