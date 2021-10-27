#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RACON } from '../../../modules/racon/main.nf' addParams( options: [:] )

workflow test_racon {
    
    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]
    assembly = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    paf      = file(params.test_data['sarscov2']['genome']['genome_paf'], checkIfExists: true)

    RACON ( input, assembly, paf)
}