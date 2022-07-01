#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTREZDIRECT_ESUMMARY } from '../../../../modules/entrezdirect/esummary/main.nf'
include { ENTREZDIRECT_XTRACT   } from '../../../../modules/entrezdirect/xtract/main.nf'

//
// Test with Assembly database
//
workflow test_entrezdirect_xtract_assembly {
    input = [
        [ id:'test_assembly' ], // meta map
        database = 'assembly',
        uid = '191021'
    ]

    pattern_in = 'DocumentSummary'
    element_in = 'SpeciesName BioprojectAccn FtpPath_GenBank'

    ENTREZDIRECT_ESUMMARY ( input )
    ENTREZDIRECT_XTRACT ( ENTREZDIRECT_ESUMMARY.out.xml_esummary, pattern_in, element_in )
}

//
// Test with Genome database
//
workflow test_entrezdirect_xtract_genome {
    input = [
        [ id:'test_genome' ], // meta map
        database = 'genome',
        uid = '768'
    ]

    pattern_in = 'DocumentSummary'
    element_in = 'TaxId Organism_Name Project_Accession Assembly_Accession'

    ENTREZDIRECT_ESUMMARY ( input )
    ENTREZDIRECT_XTRACT ( ENTREZDIRECT_ESUMMARY.out.xml_esummary, pattern_in, element_in )
}
