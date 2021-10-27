#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CMSEQ_POLYMUT as CMSEQ_POLYMUT_TEST1 ; CMSEQ_POLYMUT as CMSEQ_POLYMUT_TEST2 } from '../../../../modules/cmseq/polymut/main.nf' addParams( options: [:] )

workflow test_cmseq_polymut {
    
    input_1 = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true),
              [] ]
    

    input_2 = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true), ]

    CMSEQ_POLYMUT_TEST1 ( input_1 )
    CMSEQ_POLYMUT_TEST2 ( input_2 )
}
