#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAKECMAP_MAKECMAP } from '../../../../modules/makecmap/makecmap/main.nf'

workflow test_makecmap_makecmap {
    
    input = [
        [ id:'test'], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    MAKECMAP_MAKECMAP ( input,'bspq1' )
}
