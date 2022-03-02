#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { YARA_INDEX  } from '../../../../modules/yara/index/main.nf'
include { YARA_MAPPER } from '../../../../modules/yara/mapper/main.nf'

workflow test_yara_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    YARA_INDEX ( fasta )
    YARA_MAPPER ( input, YARA_INDEX.out.index )
}

workflow test_yara_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    YARA_INDEX ( fasta )
    YARA_MAPPER ( input, YARA_INDEX.out.index )
}
