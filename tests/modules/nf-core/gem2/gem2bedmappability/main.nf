#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GEM2_GEM2BEDMAPPABILITY   } from '../../../../../modules/nf-core/gem2/gem2bedmappability/main.nf'
include { GEM2_GEMMAPPABILITY       } from '../../../../../modules/nf-core/gem2/gemmappability/main.nf'
include { GEM2_GEMINDEXER           } from '../../../../../modules/nf-core/gem2/gemindexer/main.nf'

workflow test_gem2_gem2bedmappability {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    GEM2_GEMINDEXER ( input )

    GEM2_GEMMAPPABILITY ( GEM2_GEMINDEXER.out.index, "100" )

    GEM2_GEM2BEDMAPPABILITY ( GEM2_GEMMAPPABILITY.out.map, GEM2_GEMINDEXER.out.index )
}
