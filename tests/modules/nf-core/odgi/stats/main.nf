#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ODGI_STATS } from '../../../../../modules/nf-core/odgi/stats/main.nf'

workflow test_odgi_stats {

    input = [ [ id:'test' ],
              [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_smoothxg_gfa'], checkIfExists: true) ],
            ]

    ODGI_STATS ( input )
}
