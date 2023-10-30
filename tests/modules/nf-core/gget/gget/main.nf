#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GGET_GGET as GGET_REF } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_REF_DOWNLOAD } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_SEARCH } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_SEARCH_CSV } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_INFO } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_SEQ } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_BLAT } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_MUSCLE } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_ENRICHR } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_ARCHS4 } from '../../../../../modules/nf-core/gget/gget/main.nf'
include { GGET_GGET as GGET_PDB } from '../../../../../modules/nf-core/gget/gget/main.nf'

workflow test_gget_ref {
    
    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_REF ( input )
}

workflow test_gget_ref_download {
    
    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_REF ( input )
}

workflow test_gget_search {
    
    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_SEARCH( input )
}

workflow test_gget_search_csv {
    
    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_SEARCH_CSV ( input )
}

workflow test_gget_info {
    
    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_INFO( input )
}

workflow test_gget_seq {
    
    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_SEQ( input )
}


workflow test_gget_blat {
    
    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_BLAT( input )
}

workflow test_gget_muscle {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true)
    ]

    GGET_MUSCLE ( input )
}

workflow test_gget_enrichr {

    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_ENRICHR( input )
}

workflow test_gget_archs4 {
    
    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_ARCHS4( input )
}

workflow test_gget_pdb {
    
    input = [
        [ id:'test' ], // meta map
        []
    ]

    GGET_PDB( input )
}