#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREEBAYES } from '../../../modules/freebayes/main.nf' addParams( options: [:] )

workflow test_freebayes {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]
    reference = [file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
                 file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)]

    FREEBAYES ( input, reference )
}
