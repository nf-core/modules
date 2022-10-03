#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FCS_FCSADAPTOR } from '../../../../../modules/nf-core/fcs/fcsadaptor/main.nf'

workflow test_fcs_fcsadaptor {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    FCS_FCSADAPTOR ( input )
}
