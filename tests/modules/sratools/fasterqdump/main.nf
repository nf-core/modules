#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../../modules/untar/main.nf'
include { SRATOOLS_FASTERQDUMP } from '../../../../modules/sratools/fasterqdump/main.nf'

workflow test_sratools_fasterqdump_single_end {

    archive = [ [], file(params.test_data['sarscov2']['illumina']['SRR13255544_tar_gz'], checkIfExists: true) ]
    UNTAR ( archive )

    def input = Channel.of([ id:'test_single_end', single_end:true ])
        .combine(UNTAR.out.untar.map{ it[1] })

    SRATOOLS_FASTERQDUMP ( input )
}

workflow test_sratools_fasterqdump_paired_end {

    archive = [ [], file(params.test_data['sarscov2']['illumina']['SRR11140744_tar_gz'], checkIfExists: true) ]
    UNTAR ( archive )

    def input = Channel.of([ id:'test_paired_end', single_end:false ])
        .combine(UNTAR.out.untar.map{ it[1] })

    SRATOOLS_FASTERQDUMP ( input )
}
