#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GRAPHMAP2_INDEX } from '../../../../modules/graphmap2/index/main.nf'
include { GRAPHMAP2_ALIGN } from '../../../../modules/graphmap2/align/main.nf'

workflow test_graphmap2_align {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GRAPHMAP2_INDEX ( fasta )
    GRAPHMAP2_ALIGN ( input, fasta, GRAPHMAP2_INDEX.out.index )
}
