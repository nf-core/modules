include { KRAKEN2_ADD   } from '../../../modules/nf-core/kraken2/add/main'
include { KRAKEN2_BUILD } from '../../../modules/nf-core/kraken2/build/main'

workflow FASTA_BUILD_ADD_KRAKEN2 {

    take:
    ch_fasta              // channel: [ val(meta), fasta ]
    ch_taxonomy_names     // channel: [ names.dmp ]
    ch_taxonomy_nodes     // channel: [ nodes.dmp ]
    ch_accession2taxid    // channel: [ acc2taxidfile ]
    val_cleanintermediate // value: [ true | false ]

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

    KRAKEN2_BUILD ( KRAKEN2_BUILD.out.db, val_cleanintermediate )
    ch_versions = ch_versions.mix(KRAKEN2_BUILD.out.versions.first())

    emit:
    db = KRAKEN2_BUILD.out.db // channel: [ val(meta), [ db ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

