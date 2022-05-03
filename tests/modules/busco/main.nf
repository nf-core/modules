#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO } from '../../../modules/busco/main.nf'

// This tests genome decompression, empty input channels and data download
workflow test_busco {

    input = [
        [ id:'test', single_end:false ], // meta map
        file( params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [], // No config
    )

}

