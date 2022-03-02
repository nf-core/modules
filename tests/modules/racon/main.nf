#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RACON } from '../../../modules/racon/main.nf'

workflow test_racon {
    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['bacteroides_fragilis']['nanopore']['test_fastq_gz'], checkIfExists: true),
              file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true),
              file(params.test_data['bacteroides_fragilis']['genome']['genome_paf'], checkIfExists: true)
            ]

    RACON ( input )
}