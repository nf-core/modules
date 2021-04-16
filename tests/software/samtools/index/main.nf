#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BAI } from '../../../../software/samtools/index/main.nf' addParams( options: [:] )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CSI } from '../../../../software/samtools/index/main.nf' addParams( options: [args:'-c'] )

workflow test_samtools_index_bai {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    SAMTOOLS_INDEX_BAI ( input )
}

workflow test_samtools_index_csi {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    SAMTOOLS_INDEX_CSI ( input )
}
