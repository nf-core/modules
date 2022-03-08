#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR           } from '../../../../modules/untar/main.nf'
include { KRAKEN2_KRAKEN2 } from '../../../../modules/kraken2/kraken2/main.nf'
include { BRACKEN_BRACKEN } from '../../../../modules/bracken/bracken/main.nf'

workflow test_kraken2_bracken_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true)

    UNTAR ( db )
    KRAKEN2_KRAKEN2 ( input, UNTAR.out.untar )
    def tax_level = Channel.value('S')
    BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.txt.combine(UNTAR.out.untar).combine(tax_level) )
}

workflow test_kraken2_bracken_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2_tar_gz'], checkIfExists: true)

    UNTAR ( db )
    KRAKEN2_KRAKEN2 ( input, UNTAR.out.untar )
    def tax_level = Channel.value('S')
    BRACKEN_BRACKEN ( KRAKEN2_KRAKEN2.out.txt.combine(UNTAR.out.untar).combine(tax_level) )
}
