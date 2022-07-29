#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC            } from '../../../modules/fastqc/main.nf'
include { FASTQC as FASTQC2 } from '../../../modules/fastqc/main.nf'
include { MULTIQC           } from '../../../modules/multiqc/main.nf'

workflow test_multiqc {
    input = [
        [ id: 'test', single_end: false ],
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)]
    ]

    FASTQC  ( input )
    MULTIQC ( FASTQC.out.zip.collect { it[1] }, [[],[]] )
}

workflow test_multiqc_fn_collision {
    fqc_input = [
        [ id: 'test', single_end: false ],
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)]
    ]
    mqc_input = Channel.empty()

    FASTQC  ( fqc_input  )
    mqc_input = mqc_input.mix(FASTQC.out.zip.collect { it[1] })

    FASTQC2  ( fqc_input )
    mqc_input = mqc_input.mix(FASTQC2.out.zip.collect { it[1] })

    MULTIQC ( mqc_input, [[],[]] )
}
