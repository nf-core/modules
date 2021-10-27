#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CMSEQ_POLYMUT } from '../../../../modules/cmseq/polymut/main.nf' addParams( options: [:] )

workflow test_cmseq_polymut {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true),
              [] ]
    CMSEQ_POLYMUT ( input )
}
