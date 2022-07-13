#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GOAT_TAXONSEARCH } from '../../../../modules/goat/taxonsearch/main.nf'

//
// Test with genus name (Canis)
//
workflow test_goat_taxonsearch_genus_name {

    input = [
      [ id:'test_genus_name' ], // meta map
      taxon = 'Canis',
      []
  ]
    GOAT_TAXONSEARCH ( input )
}

//
// Test with genus (Drosophila, fruit flies) using NCBI taxonomy ID
//
workflow test_goat_taxonsearch_genus_id {

    input = [
      [ id:'test_genus_id' ], // meta map
      taxon = '7215',
      []
  ]
    GOAT_TAXONSEARCH ( input )
}
