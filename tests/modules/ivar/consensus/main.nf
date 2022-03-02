#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.save_mpileup = true
include { IVAR_CONSENSUS } from '../../../../modules/ivar/consensus/main.nf'

workflow test_ivar_consensus {
    input = [ [ id:'test'],
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    IVAR_CONSENSUS ( input, fasta )
}
