include { KRAKEN2_ADD   } from '../../../modules/nf-core/kraken2/add/main'
include { KRAKEN2_BUILD } from '../../../modules/nf-core/kraken2/build/main'
include { BRACKEN_BUILD } from '../../../modules/nf-core/bracken/build/main'

workflow FASTA_BUILD_ADD_KRAKEN2_BRACKEN {
    take:
    ch_fasta               // channel: [ val(meta), [ fasta1, fasta2, fasta3] ]
    ch_taxonomy_names      // channel: [ names.dmp ]
    ch_taxonomy_nodes      // channel: [ nodes.dmp ]
    ch_accession2taxid     // channel: [ acc2taxidfile ]
    val_cleanintermediates // value:   [ true | false ]
    ch_custom_seqid2taxid  // channel: [ seqid2taxidfile ]
    val_runbrackenbuild    // value:   [ true | false ]

    main:

    if (val_cleanintermediates && val_runbrackenbuild) {
        error("Cannot perform Kraken2 cleanup and build Bracken database. Bracken requires intermediate files")
    }
    val_cleanup = [val_cleanintermediates && !val_runbrackenbuild].any() ? true : false

    ch_versions = Channel.empty()

    KRAKEN2_ADD(ch_fasta, ch_taxonomy_names, ch_taxonomy_nodes, ch_accession2taxid, ch_custom_seqid2taxid)
    ch_versions = ch_versions.mix(KRAKEN2_ADD.out.versions.first())

    KRAKEN2_BUILD(KRAKEN2_ADD.out.db, val_cleanup)
    ch_versions = ch_versions.mix(KRAKEN2_BUILD.out.versions.first())

    if (val_runbrackenbuild) {
        BRACKEN_BUILD(KRAKEN2_BUILD.out.db)
        ch_final_db = BRACKEN_BUILD.out.db
        ch_versions = ch_versions.mix(BRACKEN_BUILD.out.versions.first())
    }
    else {
        ch_final_db = KRAKEN2_BUILD.out.db
        ch_versions = ch_versions.mix(KRAKEN2_BUILD.out.versions.first())
    }

    emit:
    db       = ch_final_db // channel: [ val(meta), [ db ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
