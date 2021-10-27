#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMCMP } from '../../../modules/bamcmp/main.nf' addParams( options: [:] )

workflow test_bamcmp {
    
    input = [ [ id:'test'],
    file('/mnt/panfs1/scratch/wsspaces/kmurat-nfcore-0/DSL2/GIT/DSL1/methylation_pipeline/results/contamination_merged/M023_CDX13_S220_C150_M1R_Input_contamination.bam', checkIfExists: true),
    file('/mnt/panfs1/scratch/wsspaces/kmurat-nfcore-0/DSL2/GIT/DSL1/methylation_pipeline/results/bam_merged/M023_CDX13_S220_C150_M1R_Input_aligned.bam', checkIfExists: true)
    
             // file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
             // file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
            ]

    BAMCMP ( input, 'bam' )
}
