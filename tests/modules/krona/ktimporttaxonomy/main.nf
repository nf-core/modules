#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONA_KTIMPORTTAXONOMY_READS }  from '../../../../modules/krona/ktimporttaxonomy/main.nf'

include { KRONA_KTIMPORTTAXONOMY_REPORT } from '../../../../modules/krona/ktimporttaxonomy/main.nf'

workflow test_krona_ktimporttaxonomy_reads {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['metagenome']['test_1.kraken2.reads.txt'], checkIfExists: true)
    ]
    taxonomy = file(params.test_data['generic']['txt']['hello'], checkIfExists: true)

    KRONA_KTIMPORTTAXONOMY_READS ( input, taxonomy )
}

workflow test_krona_ktimporttaxonomy_report {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['metagenome']['test_1.kraken2.report.txt'], checkIfExists: true)
    ]
    taxonomy = file(params.test_data['generic']['txt']['hello'], checkIfExists: true)

    KRONA_KTIMPORTTAXONOMY_REPORT ( input, taxonomy )
}
