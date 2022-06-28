#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTREZDIRECT_ESEARCH } from '../../../../modules/entrezdirect/esearch/main.nf'

//
// Test with PubMed database, using date range and spell check
//
workflow test_entrezdirect_esearch_pubmed {
    input = [
        [ id:'test_pubmed' ], // meta map
        database = 'pubmed',
        term = 'GABA+receptor'
    ]
    sort_by     =  'pub+date'
    date_type   =  'pdat'
    min_date    =  '2021/06/20'
    max_date    =  '2022/06/20'
    spell_check =  true

    ENTREZDIRECT_ESEARCH ( input, sort_by, date_type, min_date, max_date, spell_check )
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
    sort_by     =  ''
    date_type   =  ''
    min_date    =  ''
    max_date    =  ''
    spell_check =  false
    ENTREZDIRECT_ESEARCH ( input, sort_by, date_type, min_date, max_date, spell_check )
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
    sort_by     =  ''
    date_type   =  ''
    min_date    =  ''
    max_date    =  ''
    spell_check =  false
    ENTREZDIRECT_ESEARCH ( input, sort_by, date_type, min_date, max_date, spell_check )
}
