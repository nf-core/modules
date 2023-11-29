#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTREZDIRECT_ESUMMARY } from '../../../../../modules/nf-core/entrezdirect/esummary/main.nf'
include { ENTREZDIRECT_XTRACT   } from '../../../../../modules/nf-core/entrezdirect/xtract/main.nf'

//
// Test with Assembly database
//
workflow test_entrezdirect_xtract_assembly {

    input = [
        [ id:'test_assembly' ], // meta map
        uid = '191021',
        []
    ]
    database = 'assembly'

    pattern_in = 'DocumentSummary'
    element_in = 'SpeciesName BioprojectAccn FtpPath_GenBank'
    delim = ","

    ENTREZDIRECT_ESUMMARY ( input, database )
    ENTREZDIRECT_XTRACT ( ENTREZDIRECT_ESUMMARY.out.xml, pattern_in, element_in, delim )
}

//
// Test with Genome database
//
workflow test_entrezdirect_xtract_genome {

    input = [
        [ id:'test_genome' ], // meta map
        uid = '768',
        []
    ]
    database = 'genome'

    pattern_in = 'DocumentSummary'
    element_in = 'TaxId Organism_Name Project_Accession Assembly_Accession'
    delim = ","

    ENTREZDIRECT_ESUMMARY ( input, database )
    ENTREZDIRECT_XTRACT ( ENTREZDIRECT_ESUMMARY.out.xml, pattern_in, element_in, delim )
}
