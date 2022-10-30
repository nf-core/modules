#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONA_KTIMPORTTAXONOMY as KRONA_KTIMPORTTAXONOMY_READS  } from '../../../../../modules/nf-core/krona/ktimporttaxonomy/main.nf'
include { KRONA_KTIMPORTTAXONOMY as KRONA_KTIMPORTTAXONOMY_REPORT } from '../../../../../modules/nf-core/krona/ktimporttaxonomy/main.nf'

workflow test_krona_ktimporttaxonomy_reads {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['metagenome']['classified_reads_assignment'], checkIfExists: true)
    ]
    taxonomy = file(params.test_data['sarscov2']['metagenome']['krona_taxonomy'], checkIfExists: true)

    KRONA_KTIMPORTTAXONOMY_READS ( input, taxonomy )
}

workflow test_krona_ktimporttaxonomy_report {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['metagenome']['kraken_report'], checkIfExists: true)
    ]
    taxonomy = file(params.test_data['sarscov2']['metagenome']['krona_taxonomy'], checkIfExists: true)

    KRONA_KTIMPORTTAXONOMY_REPORT ( input, taxonomy )
}
