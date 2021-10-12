#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../../modules/untar/main.nf'
include { SRATOOLS_FASTERQDUMP } from '../../../../modules/sratools/fasterqdump/main.nf' addParams( options: [:] )

workflow test_sratools_fasterqdump_single_end {

    archive = file(params.test_data['sarscov2']['illumina']['test_single_end_sra_tar_gz'], checkIfExists: true)
    UNTAR ( archive )

    def input = Channel.of([ id:'test_single_end', single_end:true ])
        .combine(UNTAR.out.untar)

    SRATOOLS_FASTERQDUMP ( input )
}

workflow test_sratools_fasterqdump_paired_end {

    archive = file(params.test_data['sarscov2']['illumina']['test_paired_end_sra_tar_gz'], checkIfExists: true)
    UNTAR ( archive )

    def input = Channel.of([ id:'test_paired_end', single_end:false ])
        .combine(UNTAR.out.untar)

    SRATOOLS_FASTERQDUMP ( input )
}
