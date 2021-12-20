#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BAI  } from '../../../../modules/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CRAI } from '../../../../modules/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CSI  } from '../../../../modules/samtools/index/main.nf'

workflow test_samtools_index_bai {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    SAMTOOLS_INDEX_BAI ( input )
}

workflow test_samtools_index_crai {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_cram'], checkIfExists: true)
            ]

    SAMTOOLS_INDEX_CRAI ( input )
}

workflow test_samtools_index_csi {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    SAMTOOLS_INDEX_CSI ( input )
}
