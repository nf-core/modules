include { MMSEQS_CREATEDB  } from '../../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_DATABASES } from '../../../modules/nf-core/mmseqs/databases/main'
include { MMSEQS_TAXONOMY  } from '../../../modules/nf-core/mmseqs/taxonomy/main'
include { MMSEQS_CREATETSV } from '../../../modules/nf-core/mmseqs/createtsv/main'

workflow MMSEQS_CONTIG_TAXONOMY {

    take:
    contigs            // channel: tuple val(meta), path(contigs)
    mmseqs_databases   // channel: path(mmseqs2 local db)
    databases_id       // channel: [mmseqs2_db_id]

    main:

    ch_versions               = Channel.empty()
    ch_mmseqs_db              = Channel.empty()
    ch_taxonomy_querydb       = Channel.empty()
    ch_taxonomy_querydb_taxdb = Channel.empty()
    ch_taxonomy_tsv           = Channel.empty()

    // Download the ref db if not supplied by user
    // MMSEQS_DATABASE
    if ( !mmseqs_databases.empty ) {
        ch_mmseqs_db = Channel
            .fromPath( mmseqs_databases )
            .first()
    } else {
        MMSEQS_DATABASES ( databases_id )
        ch_versions  = ch_versions.mix( MMSEQS_DATABASES.out.versions )
        ch_mmseqs_db = ( MMSEQS_DATABASES.out.database )
    }

    // Create db for query contigs, assign taxonomy and convert to table format
    // MMSEQS_CREATEDB
    MMSEQS_CREATEDB ( contigs )
    ch_versions         = ch_versions.mix( MMSEQS_CREATEDB.out.versions )
    ch_taxonomy_querydb = MMSEQS_CREATEDB.out.db

    // MMSEQS_TAXONOMY
    MMSEQS_TAXONOMY ( ch_taxonomy_querydb, ch_mmseqs_db )
    ch_versions               = ch_versions.mix( MMSEQS_TAXONOMY.out.versions )
    ch_taxonomy_querydb_taxdb = MMSEQS_TAXONOMY.out.db_taxonomy

    // MMSEQS_CREATETSV
    MMSEQS_CREATETSV ( ch_taxonomy_querydb_taxdb, [[:],[]], ch_taxonomy_querydb )
    ch_versions     = ch_versions.mix( MMSEQS_CREATETSV.out.versions )
    ch_taxonomy_tsv = MMSEQS_CREATETSV.out.tsv

    emit:
    taxonomy    = ch_taxonomy_tsv           // channel: [ val(meta), tsv ]
    db_mmseqs   = ch_mmseqs_db              // channel: [ val(meta), mmseqs_database ]
    db_taxonomy = ch_taxonomy_querydb_taxdb // channel: [ val(meta), db_taxonomy ]
    db_contig   = ch_taxonomy_querydb       // channel: [ val(meta), db ]
    versions    = ch_versions               // channel: [ versions.yml ]
}
