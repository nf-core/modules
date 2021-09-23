#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEACR_CALLPEAKNOCTRL } from '../../../../modules/seacr/callpeaknoctrl/main.nf' addParams( options: [ args:'norm stringent' ] )

workflow test_seacr_callpeaknoctrl {
    input = [ [ id:'test_1'],
              file(params.test_data['homo_sapiens']['illumina']['cutandrun_bedgraph_test_1'], checkIfExists: true)
            ]

    SEACR_CALLPEAKNOCTRL ( input, 0.05 )
}
