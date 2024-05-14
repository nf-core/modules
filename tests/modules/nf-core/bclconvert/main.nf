#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCLCONVERT } from '../../../../modules/nf-core/bclconvert/main.nf'

workflow test_bclconvert {
    //TODO use new test dataset when available, see https://github.com/nf-core/test-datasets/issues/996
    ch_flowcell = Channel.value([
            [id:'test', lane:1 ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell_samplesheet'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell'], checkIfExists: true)])

    ch_flowcell
    .multiMap { meta, ss, run ->
        samplesheet: [meta, ss]
        tar: [meta, run]
    }.set{ ch_fc_split }

    ch_flowcell_merge = ch_fc_split.samplesheet.join( ch_fc_split.tar )

    BCLCONVERT (ch_flowcell_merge)
}
