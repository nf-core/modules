#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.save_mpileup = true
include { IVAR_CONSENSUS } from '../../../../software/ivar/consensus/main.nf' addParams( [ options: [args2: '-aa -A -d 0 -Q 0'] ] )

workflow test_ivar_consensus {
    input = [ [ id:'test'],
                file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    
    IVAR_CONSENSUS ( input, fasta )
}
