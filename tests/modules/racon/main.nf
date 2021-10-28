#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RACON } from '../../../modules/racon/main.nf' addParams( options: [:] )

workflow test_racon {
    input = [ [ id:'test', single_end:true ], // meta map
             [ file(params.test_data['bacteroides_fragilis']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]
    assembly = file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    paf      = file(params.test_data['bacteroides_fragilis']['genome']['genome_paf'], checkIfExists: true)

    RACON ( input, assembly, paf)
}