#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ODGI_SQUEEZE } from '../../../../../modules/nf-core/odgi/squeeze/main.nf'

workflow test_odgi_squeeze {
    og_file = file(params.test_data['homo_sapiens']['pangenome']['odgi']['pangenome_og'], checkIfExists: true)
    gfa_file = file(params.test_data['homo_sapiens']['pangenome']['pangenome_smoothxg_gfa'], checkIfExists: true)
    
    ch_input = Channel.from([ 
        [ id:'test', single_end:false ], // meta map
        [ og_file, gfa_file ]
    ])

    ODGI_SQUEEZE ( ch_input )
}
