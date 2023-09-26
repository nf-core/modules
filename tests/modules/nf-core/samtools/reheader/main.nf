#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_REHEADER as SAMTOOLS_REHEADER_RGDEL} from '../../../../../modules/nf-core/samtools/reheader/main.nf'
include { SAMTOOLS_REHEADER as SAMTOOLS_REHEADER_CHRDEL} from '../../../../../modules/nf-core/samtools/reheader/main.nf'

workflow test_samtools_reheader_chrdel {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    SAMTOOLS_REHEADER_CHRDEL ( input )
}

workflow test_samtools_reheader_rgdel {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    SAMTOOLS_REHEADER_RGDEL ( input )
}