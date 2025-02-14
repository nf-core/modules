#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LEEHOM } from '../../../../modules/nf-core/leehom/main.nf'
include { SAMTOOLS_VIEW } from '../../../../modules/nf-core/samtools/view/main.nf'

workflow test_leehom_bam {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            []
        ]

    SAMTOOLS_VIEW ( input, [] , [], [])
    LEEHOM ( SAMTOOLS_VIEW.out.bam )
}

workflow test_leehom_se_fq {

    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]

    LEEHOM ( input )
}

workflow test_leehom_pe_fq {

    input = [ [ id:'test', single_end:false ], // meta map
              [
                  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
              ] ]

    LEEHOM ( input )
}
