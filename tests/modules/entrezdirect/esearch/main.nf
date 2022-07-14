#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTREZDIRECT_ESEARCH as ENTREZDIRECT_ESEARCHP  } from '../../../../modules/entrezdirect/esearch/main.nf'
include { ENTREZDIRECT_ESEARCH                           } from '../../../../modules/entrezdirect/esearch/main.nf'

//
// Test with PubMed database, using date range and spell check,
// see nextflow.config file for optional definition(ext.args)
//
workflow test_entrezdirect_esearch_pubmed {
    input = [
        [ id:'test_pubmed' ], // meta map
        database = 'pubmed',
        term = 'selective serotonin reuptake inhibitor'
    ]

    ENTREZDIRECT_ESEARCHP ( input )
}

//
// Test with Genome database and species; no date range, no spell check
//
workflow test_entrezdirect_esearch_genome {
    input = [
        [ id:'test_genome' ], // meta map
        database = 'genome',
        term = 'Danio+rerio'
    ]

    ENTREZDIRECT_ESEARCH ( input )
}

//
// Test with Assembly database and GenBank accession; no date range, no spell check
//
workflow test_entrezdirect_esearch_assembly {
    input = [
        [ id:'test_assembly' ], // meta map
        database = 'assembly',
        term = 'GCA_000001635.9'
    ]

    ENTREZDIRECT_ESEARCH ( input )
}
