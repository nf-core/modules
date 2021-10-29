#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MEDAKA } from '../../../modules/medaka/main.nf' addParams( options: [:] )

workflow test_medaka {
    
    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['bacteroides_fragilis']['nanopore']['test_fastq_gz'], checkIfExists: true),
              file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
            ]
    
    MEDAKA ( input )
}
