#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMCMP } from '../../../modules/bamcmp/main.nf' 

workflow test_bamcmp {
    
    input = [
    file('/mnt/panfs1/scratch/wsspaces/kmurat-nfcore-0/DSL2/GIT/DSL1/methylation_pipeline/results/contamination_merged/M023_CDX13_S220_C150_M1R_Input_contamination.bam', checkIfExists: true),
    file('/mnt/panfs1/scratch/wsspaces/kmurat-nfcore-0/DSL2/GIT/DSL1/methylation_pipeline/results/bam_merged/M023_CDX13_S220_C150_M1R_Input_aligned.bam.bam', checkIfExists: true)
            ]

    BAMCMP ( input )
}
