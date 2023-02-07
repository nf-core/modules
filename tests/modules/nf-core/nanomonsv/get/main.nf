#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOMONSV_GET } from '../../../../../modules/nf-core/nanomonsv/get/main.nf'

workflow test_nanomonsv_get {

    def data_path = '/Users/arthurgymer/workbench/nanomonsv/tests' //'https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/nanomonsv'
    input = channel.of([
        [ id:'test' ], // meta map
        file("${data_path}/resource/bam/test_tumor.bam", checkIfExists: true),
        file("${data_path}/resource/bam/test_tumor.bam.bai", checkIfExists: true),
        file("${data_path}/data/test_tumor/*.gz*"),
        file("${data_path}/resource/bam/test_ctrl.bam", checkIfExists: true),
        file("${data_path}/resource/bam/test_ctrl.bam.bai", checkIfExists: true),
        file("${data_path}/data/test_ctrl/*.gz*"),
    ])
    ref = file('s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa')

    NANOMONSV_GET ( input, ref )
}
