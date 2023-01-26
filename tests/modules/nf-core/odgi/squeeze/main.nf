#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ODGI_SQUEEZE } from '../../../../../modules/nf-core/odgi/squeeze/main.nf'

workflow test_odgi_squeeze {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['pangenome']['odgi']['pangenome_og'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['pangenome']['pangenome_smoothxg_gfa'], checkIfExists: true)]
            ]

    ODGI_SQUEEZE ( input )
}
