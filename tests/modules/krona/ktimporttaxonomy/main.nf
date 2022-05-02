#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONA_KTIMPORTTAXONOMY as TAXONOMY_READS  } from '../../../../modules/krona/ktimporttaxonomy/main.nf'

include { KRONA_KTIMPORTTAXONOMY as TAXONOMY_REPORT } from '../../../../modules/krona/ktimporttaxonomy/main.nf'

include { KRONA_KRONADB          as DOWNLOAD_DB     } from '../../../../modules/krona/kronadb/main.nf'

workflow test_krona_ktimporttaxonomy_reads {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['metagenome']['classified_reads_assignment'], checkIfExists: true)
    ]
    // taxonomy = file(params.test_data['generic']['txt']['hello'], checkIfExists: true)

    DOWNLOAD_DB()
    TAXONOMY_READS ( input, DOWNLOAD_DB.out.db )
}

workflow test_krona_ktimporttaxonomy_report {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['metagenome']['report'], checkIfExists: true)
    ]
    taxonomy = file(params.test_data['generic']['txt']['hello'], checkIfExists: true)

    DOWNLOAD_DB()
    TAXONOMY_REPORT ( input, DOWNLOAD_DB.out.db )
}
