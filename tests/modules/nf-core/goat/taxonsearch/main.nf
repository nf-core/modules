#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GOAT_TAXONSEARCH } from '../../../../../modules/nf-core/goat/taxonsearch/main.nf'

//
// Test with genus name (Canis)
//
workflow test_goat_taxonsearch_single_species {

    input = [
      [ id:'test_single_species' ], // meta map
      taxon = 'Meles meles',
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

//
// Test with multiple species from a taxa file
//
workflow test_goat_taxonsearch_species {

    input = [
      [ id:'test_species' ], // meta map
      taxon = '',
      file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/goat_taxonsearch/taxonomy_ids.txt', checkIfExists: true)
  ]
    GOAT_TAXONSEARCH ( input )
}
