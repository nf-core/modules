#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../../modules/nf-core/bowtie2/build/main.nf'
include { RENAME_PATH } from '../../../../../modules/nf-core/fastqscreen/buildfromindex/main.nf'
include { FASTQSCREEN_BUILDFROMINDEX } from '../../../../../modules/nf-core/fastqscreen/buildfromindex/main.nf'

workflow test_fastqscreen_buildfromindex {

   // fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    // fasta = [
    //     [ id:'test'],
    //     file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    // ]



    // input = [
    //     [ id:'test', single_end:false ], // meta map
    //     file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    // ]



    ch_fasta = Channel.from([
        [[id: 'test1'], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)],
        [[id: 'test21'], file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)]

        ])


    BOWTIE2_BUILD ( ch_fasta )

    RENAME_PATH(BOWTIE2_BUILD.out.index)

    RENAME_PATH.out.restaged.view()

    FASTQSCREEN_BUILDFROMINDEX (RENAME_PATH.out.restaged.map{it -> it[0]}.collect(),
                                             RENAME_PATH.out.restaged.map{it -> it[1]}.collect())

    FASTQSCREEN_BUILDFROMINDEX.out.database.view()
    FASTQSCREEN_BUILDFROMINDEX.out.versions.view()
}
