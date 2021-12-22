#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRA_FASTQ_SRATOOLS } from '../../../../subworkflows/nf-core/sra_fastq_sratools/main.nf'

workflow test_sra_fastq_sratools_single_end {
    input = [
        [ id:'test_single_end', single_end:true ], // meta map
        'SRR13255544'
    ]

    SRA_FASTQ_SRATOOLS ( input )
}

workflow test_sra_fastq_sratools_paired_end {
    input = [
        [ id:'test_paired_end', single_end:false ], // meta map
        'SRR11140744'
    ]

    SRA_FASTQ_SRATOOLS ( input )
}
