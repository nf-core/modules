#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAPDAMAGE2 } from '../../../modules/mapdamage2/main.nf'

workflow test_mapdamage2 {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) 
            ]                            
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    
    MAPDAMAGE2 ( input, fasta )
}
