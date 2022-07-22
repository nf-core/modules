#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATLAS_RECAL } from '../../../../modules/atlas/recal/main.nf'

workflow test_atlas_recal {
    
    input = [
        [ id:'test', single_end:false ], // meta map
    //file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    file("/mnt/archgen/DAG_GL/s_m/PYL004_mergedReads.bam", checkIfExists:true), 
    file("/mnt/archgen/DAG_GL/s_m/PYL004_mergedReads.bam.bai", checkIfExists:true), 
    file("/mnt/archgen/DAG_GL/infos/PYL004_PMD_input_Empiric.txt", checkIfExists:true)
    ]

    input2 = [
    file("/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa", checkIfExists:true), 
    file("/mnt/archgen/users/gnecchi/Jena/HISTOGENES/analysis/88_mammals.epo_low_coverage.10M_GRCh37.masked.bed", checkIfExists:true)
    ]

    ATLAS_RECAL ( input, input2 )
}


