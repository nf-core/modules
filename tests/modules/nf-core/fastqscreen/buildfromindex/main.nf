#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../../modules/nf-core/bowtie2/build/main.nf'
include { FASTQSCREEN_BUILDFROMINDEX } from '../../../../../modules/nf-core/fastqscreen/buildfromindex/main.nf'

workflow test_fastqscreen_buildfromindex {


    ch_fasta = Channel.from([
        [[id: 'sarscov2'], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)],
        [[id: 'homo_sapiens21'], file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)]
        ])


    BOWTIE2_BUILD ( ch_fasta )

    FASTQSCREEN_BUILDFROMINDEX (BOWTIE2_BUILD.out.index.map{it -> it[0]}.collect(),
                                BOWTIE2_BUILD.out.index.map{it -> it[1]}.collect())

    FASTQSCREEN_BUILDFROMINDEX.out.database.view()
    FASTQSCREEN_BUILDFROMINDEX.out.versions.view()
}
