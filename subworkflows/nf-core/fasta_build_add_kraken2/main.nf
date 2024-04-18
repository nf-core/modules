include { KRAKEN2_ADD   } from '../../../modules/nf-core/kraken2/add/main'
include { KRAKEN2_BUILD } from '../../../modules/nf-core/kraken2/build/main'
include { BRACKEN_BUILD } from '../../../modules/nf-core/bracken/build/main'

workflow FASTA_BUILD_ADD_KRAKEN2 {

    take:
    ch_fasta              // channel: [ val(meta), fasta ]
    ch_taxonomy_names     // channel: [ names.dmp ]
    ch_taxonomy_nodes     // channel: [ nodes.dmp ]
    ch_accession2taxid    // channel: [ acc2taxidfile ]
    val_cleanintermediate // value: [ true | false ]
    run_brackenbuild      // value: [ true | false ]

    main:

    ch_versions = Channel.empty()

    ch_fastas_for_kraken2add = ch_fasta
                                .map {
                                    meta, fasta ->

                                    [[id: 'db'], fasta]
                                }
                                .groupTuple()

    KRAKEN2_ADD ( ch_fastas_for_kraken2add, ch_taxonomy_names, ch_taxonomy_nodes, ch_accession2taxid )
    ch_versions = ch_versions.mix(KRAKEN2_ADD.out.versions.first())

    KRAKEN2_BUILD ( KRAKEN2_ADD.out.db, val_cleanintermediate )
    
    if ( run_brackenbuild ) {
        BRACKEN_BUILD ( KRAKEN2_BUILD.out.db )
        out_db = BRACKEN_BUILD.out.db
        ch_versions = ch_versions.mix(BRACKEN_BUILD.out.versions.first())
    }
    else {
        out_db = KRAKEN2_BUILD.out.db
        ch_versions = ch_versions.mix(KRAKEN2_BUILD.out.versions.first())
    }

    emit:
    db = out_db               // channel: [ val(meta), [ db ] ]
    versions = ch_versions    // channel: [ versions.yml ]
}

