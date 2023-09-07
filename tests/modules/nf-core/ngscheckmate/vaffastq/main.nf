#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_BUILD } from '../../../../../modules/nf-core/bowtie/build/main.nf'
include { BEDTOOLS_MAKEWINDOWS } from '../../../../../modules/nf-core/bedtools/makewindows/main.nf'
include { NGSCHECKMATE_PATTERNGENERATOR } from '../../../../../modules/nf-core/ngscheckmate/patterngenerator/main.nf'
include { NGSCHECKMATE_FASTQ } from '../../../../../modules/nf-core/ngscheckmate/fastq/main.nf'
include { NGSCHECKMATE_VAFFASTQ } from '../../../../../modules/nf-core/ngscheckmate/vaffastq/main.nf'

workflow test_ngscheckmate_vaffastq {

    input = [
        [
            [ id:'test1', single_end:false ], // meta map
            [
                file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
            ]
        ],
        [
            [ id:'test2', single_end:false ], // meta map
            [
                file(params.test_data['sarscov2']['illumina']['test2_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_2_fastq_gz'], checkIfExists: true)
            ]
        ],
    ]

    input_bed = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]

    fasta    = [
        [ id: 'sarscov2' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]

    bowtie_fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    BOWTIE_BUILD ( bowtie_fasta )

    NGSCHECKMATE_PATTERNGENERATOR ( input_bed, fasta, BOWTIE_BUILD.out.index )
    NGSCHECKMATE_FASTQ ( Channel.from(input), NGSCHECKMATE_PATTERNGENERATOR.out.pt.map{it[1]} )
    NGSCHECKMATE_VAFFASTQ ( NGSCHECKMATE_FASTQ.out.vaf.map{it[1]}.collect() )
}
