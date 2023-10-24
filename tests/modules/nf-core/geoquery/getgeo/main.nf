#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GEOQUERY_GETGEO } from '../../../../../modules/nf-core/geoquery/getgeo/main.nf'

workflow test_geoquery_getgeo {
    
    input = [[ id:'test' ], // meta map
        'GSE50790'
    ]

    GEOQUERY_GETGEO ( input )
}

workflow test_geoquery_getgeo_with_metacols {
    
    input = [[ id:'test' ], // meta map
        'GSE50790'
    ]

    GEOQUERY_GETGEO ( input )
}
