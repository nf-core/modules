#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CMSEQ_POLYMUT } from '../../../../modules/cmseq/polymut/main.nf'

workflow test_cmseq_polymut_1 {

    input_1 = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              [],
              file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true),
              [] ]

    CMSEQ_POLYMUT( input_1 )
    
}

workflow test_cmseq_polymut_2 {
    input_2 = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true),
              [] ]

    CMSEQ_POLYMUT( input_2 )
}

workflow test_cmseq_polymut_3 {
    input_3 = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true), ]

    CMSEQ_POLYMUT( input_3 )
}
    
