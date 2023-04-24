#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_BUILD } from '../../../../../modules/nf-core/bowtie/build/main.nf'
include { RENAME_PATH } from '../../../../../modules/nf-core/fastqscreen/buildfromindex/main.nf'
include { FASTQSCREEN_BUILDFROMINDEX } from '../../../../../modules/nf-core/fastqscreen/buildfromindex/main.nf'

workflow test_fastqscreen_buildfromindex {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BOWTIE_BUILD ( fasta )

    // input = [
    //     [ id:'test', single_end:false ], // meta map
    //     file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    // ]

    // ch_fasta = Channel.from([
    //     [[id: 'test1'], file(params.index_fasta1, checkIfExists: true)],
    //     [[id:'test2'], file(params.index_fasta2, checkIfExists: true)]
    //     ])

    RENAME_PATH(BOWTIE_BUILD.out.index)

    FASTQ_SCREEN_BUILD_DATABASE_FROM_INDEXES(RENAME_PATH.out.restaged.map{it -> it[0]['id']}.collect(),
                                             RENAME_PATH.out.restaged.map{it -> it[1]}.collect())
}
