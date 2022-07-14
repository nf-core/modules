#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCL2FASTQ } from '../../../modules/bcl2fastq/main.nf'
include { UNTAR      } from '../../../modules/untar/main.nf'

workflow test_bcl2fastq {
    ch_flowcell = Channel.value([
            [id:'test', lane:1 ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell_samplesheet'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell'], checkIfExists: true)])

    ch_flowcell
    .multiMap { meta, ss, run ->
        samplesheet: [meta, ss]
        tar: [meta, run]
    }.set{ ch_fc_split }

    ch_flowcell_untar = ch_fc_split.samplesheet.join( UNTAR ( ch_fc_split.tar ).untar )

    BCL2FASTQ (ch_flowcell_untar)
}
