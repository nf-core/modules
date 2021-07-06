#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.save_mpileup = true
include { IVAR_CONSENSUS } from '../../../../software/ivar/consensus/main.nf' addParams( [ options: [args2: '-aa -A -d 0 -Q 0'] ] )

workflow test_ivar_consensus {
    input = [ [ id:'test'],
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    IVAR_CONSENSUS ( input, fasta )
}
