include { KRAKEN2_ADD   } from '../../../modules/nf-core/kraken2/add/main'
include { KRAKEN2_BUILD } from '../../../modules/nf-core/kraken2/build/main'

workflow FASTA_BUILD_ADD_KRAKEN2 {
    take:
    ch_fasta // channel: [ val(meta), [ fasta1, fasta2, fasta3] ]
    ch_taxonomy_names // channel: [ names.dmp ]
    ch_taxonomy_nodes // channel: [ nodes.dmp ]
    ch_accession2taxid // channel: [ acc2taxidfile ]
    val_cleanintermediate // value: [ true | false ]
    ch_custom_seqid2taxid // channel: [ seqid2taxidfile ]

    main:

    ch_versions = channel.empty()

    KRAKEN2_ADD(ch_fasta, ch_taxonomy_names, ch_taxonomy_nodes, ch_accession2taxid, ch_custom_seqid2taxid)
    ch_versions = ch_versions.mix(KRAKEN2_ADD.out.versions.first())

    KRAKEN2_BUILD(KRAKEN2_ADD.out.library_added_files, KRAKEN2_ADD.out.seqid2taxid_map.ifEmpty([[:], []]), KRAKEN2_ADD.out.taxonomy_files, val_cleanintermediate)
    ch_versions = ch_versions.mix(KRAKEN2_BUILD.out.versions.first())

    emit:
    db           = KRAKEN2_BUILD.out.db // channel: [ val(meta), [ db ] ]
    db_separated = KRAKEN2_BUILD.out.db_separated // channel: [ val(meta), path(k2d_files), path(map_files), path(added_files), path(taxonomy_files) ]
    versions     = ch_versions // channel: [ versions.yml ]
}
