#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RTGTOOLS_FORMAT } from '../../../../../modules/nf-core/rtgtools/format/main.nf'

workflow test_rtgtools_format_fasta {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        [],
        []
    ]

    RTGTOOLS_FORMAT ( input )
}

workflow test_rtgtools_format_single_end_fastq {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        [],
        []
    ]

    RTGTOOLS_FORMAT ( input )
}

workflow test_rtgtools_format_paired_end_fastq {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true),
        []
    ]

    RTGTOOLS_FORMAT ( input )
}

workflow test_rtgtools_format_bam {

    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true),
        []
    ])

    sam_rg = Channel.of("@RG\tID:READGROUP1\tSM:SAMPLE\tPL:ILLUMINA")
                    .collectFile(name:'sam_rg.txt')

    RTGTOOLS_FORMAT ( input.combine(sam_rg) )
}
