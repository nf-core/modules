#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { YARA_INDEX } from '../../../../software/yara/index/main.nf' addParams( options: ['args': '-e 3'] )
include { YARA_MAPPER } from '../../../../software/yara/mapper/main.nf' addParams( options: ['args': '-e 3'] )

workflow test_yara_single_end {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    input = [ [ id:'test', single_end:true ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]


    YARA_INDEX ( fasta )
    YARA_MAPPER ( input, YARA_INDEX.out.index )
}

workflow test_yara_paired_end {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ] ]

    YARA_INDEX ( fasta )
    YARA_MAPPER ( input, YARA_INDEX.out.index )
}
