#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GET_NANOLYSE_FASTA;
          NANOLYSE           } from '../../../software/nanolyse/main.nf'    addParams( options: [:] )

workflow test_nanolyse {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)]
            ]

    GET_NANOLYSE_FASTA()
    NANOLYSE ( input, GET_NANOLYSE_FASTA.out.nanolyse_fasta )
}
