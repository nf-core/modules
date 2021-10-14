#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAXBIN2 } from '../../../modules/maxbin2/main.nf' addParams( options: [:] )

workflow test_maxbin2 {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['genome_fasta'], checkIfExists: true),
              file(params.test_data['candidatus_portiera_aleyrodidarum']['illumina']['test_1_fastq_gz'], checkIfExists: true),
              [] ]

    MAXBIN2 ( input )
}
