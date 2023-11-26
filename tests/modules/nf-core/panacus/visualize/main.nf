#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PANACUS_VISUALIZE } from '../../../../../modules/nf-core/panacus/visualize/main.nf'

workflow test_panacus_visualize {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['pangenome']['pangenome_panacus_tsv'], checkIfExists: true)
    ]

    PANACUS_VISUALIZE ( input )
}
