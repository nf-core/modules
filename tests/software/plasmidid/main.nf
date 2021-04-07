#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIA     } from '../../../software/minia/main.nf'     addParams ( options: [:] )
include { PLASMIDID } from '../../../software/plasmidid/main.nf' addParams ( options: ['args' : '-k 0.8'] )

workflow test_plasmidid {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    MINIA ( input )

    PLASMIDID ( MINIA.out.contigs, fasta )
}
