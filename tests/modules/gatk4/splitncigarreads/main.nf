#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SPLITNCIGARREADS } from '../../../../modules/gatk4/splitncigarreads/main.nf'

workflow test_gatk4_splitncigarreads {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ] 
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_SPLITNCIGARREADS ( input, fasta, fai, dict )
}
