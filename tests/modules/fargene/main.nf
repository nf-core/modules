#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../modules/gunzip/main.nf'
include { FARGENE } from '../../../modules/fargene/main.nf'

workflow test_fargene {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true) ]
    hmm_model = 'class_a'

    GUNZIP ( input )
    FARGENE ( GUNZIP.out.gunzip, hmm_model )
}
