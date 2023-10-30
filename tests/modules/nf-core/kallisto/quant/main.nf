#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTO_QUANT } from '../../../../../modules/nf-core/kallisto/quant/main.nf'
include { KALLISTO_INDEX } from '../../../../../modules/nf-core/kallisto/index/main.nf'

fasta = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
gtf   = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

workflow test_kallisto_quant {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    KALLISTO_INDEX ( fasta )
    KALLISTO_QUANT (
        input,
        KALLISTO_INDEX.out.idx,
        [],
        []
    )
}

workflow test_kallisto_quant_single_end {
    input = [
        [ id:'test', single_end:true ],
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    KALLISTO_INDEX ( fasta )
    KALLISTO_QUANT (
        input,
        KALLISTO_INDEX.out.idx,
        [],
        []
    )
}
