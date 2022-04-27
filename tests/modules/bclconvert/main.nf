#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCLCONVERT } from '../../../modules/bclconvert/main.nf'
include { UNTAR      } from '../../../modules/untar/main.nf'

workflow test_bclconvert {
    ch_flowcell = Channel.from([[ id:'test', lane:1 ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell_samplesheet'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell'], checkIfExists: true)])

    ch_flowcell
    .multiMap { meta, ss, run ->
        samplesheet: [meta, ss]
        tar: [meta, tar]
    }.set{ ch_fc_split }

    ch_flowcell_untar = ch_fc_split.samplesheet.join( UNTAR ( ch_fc_split.tar ).untar )

    BCLCONVERT (ch_flowcell_untar)
}
