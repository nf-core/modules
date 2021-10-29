#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMCMP } from '../../../modules/bamcmp/main.nf' 

workflow test_bamcmp {

    input = [[ id:'test' ],
    file('/mnt/panfs1/scratch/wsspaces/kmurat-nfcore-0/DSL2/GIT/Input_aligned.bam', checkIfExists: true),
    file('/mnt/panfs1/scratch/wsspaces/kmurat-nfcore-0/DSL2/GIT/Input_contamination.bam', checkIfExists: true)]


    BAMCMP ( input )
}
