#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIASM } from '../../../modules/miniasm/main.nf'

workflow test_miniasm {
    
    input = [ [ id:'test', single_end:true ], // meta map
             file(params.test_data['bacteroides_fragilis']['nanopore']['test_fastq_gz'], checkIfExists: true),
             file(params.test_data['bacteroides_fragilis']['nanopore']['overlap_paf'], checkIfExists: true)
            ]
    
    MINIASM ( input )
}
