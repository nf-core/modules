#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRA_FASTQ } from '../../../../subworkflows/nf-core/sra_fastq/main.nf' addParams( [:] )

workflow test_sra_fastq_single_end {
    input = [
        [ id:'test_single_end', single_end:true ], // meta map
        'SRR13255544'
    ]

    SRA_FASTQ ( input )
}

workflow test_sra_fastq_paired_end {
    input = [
        [ id:'test_paired_end', single_end:false ], // meta map
        'SRR11140744'
    ]

    SRA_FASTQ ( input )
}
