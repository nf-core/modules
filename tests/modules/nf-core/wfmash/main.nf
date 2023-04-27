#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WFMASH } from '../../../../modules/nf-core/wfmash/main.nf'

workflow test_wfmash {
    gzi = file(params.test_data['homo_sapiens']['pangenome']['pangenome_fa_bgzip_gzi'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['pangenome']['pangenome_fa_bgzip_fai'], checkIfExists: true)

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_fa_bgzip'], checkIfExists: true) ],
              [ ], // empty paf input
              gzi,
              fai
            ]
    query_self = true

    WFMASH ( input, query_self, [])
}
