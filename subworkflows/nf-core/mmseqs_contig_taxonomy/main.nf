include { MMSEQS_CREATEDB  } from '../../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_DATABASES } from '../../../modules/nf-core/mmseqs/databases/main'
include { MMSEQS_TAXONOMY  } from '../../../modules/nf-core/mmseqs/taxonomy/main'
include { MMSEQS_CREATETSV } from '../../../modules/nf-core/mmseqs/createtsv/main'

workflow MMSEQS_CONTIG_TAXONOMY {

    take:
    ch_contigs                      // channel: tuple val(meta), path(contigs)
    ch_mmseqs_databases_localpath   // channel: path(mmseqs2_db)
    ch_mmseqs_databases_id          // channel: [mmseqs2_db_id]

    main:

    ch_versions               = Channel.empty()
    ch_mmseqs_db              = Channel.empty()
    ch_taxonomy_querydb       = Channel.empty()
    ch_taxonomy_querydb_taxdb = Channel.empty()
    ch_taxonomy_tsv           = Channel.empty()

    // Download the ref db if not supplied by user
    // MMSEQS_DATABASE
    if ( !ch_mmseqs_databases_localpath.empty ) {
        ch_mmseqs_db = Channel
            .fromPath( ch_mmseqs_databases_localpath )
            .first()
    } else {
        MMSEQS_DATABASES ( ch_mmseqs_databases_id )
        ch_versions  = ch_versions.mix( MMSEQS_DATABASES.out.versions )
        ch_mmseqs_db = ( MMSEQS_DATABASES.out.database )
    }

    // Create db for query contigs, assign taxonomy and convert to table format
    // MMSEQS_CREATEDB
    MMSEQS_CREATEDB ( ch_contigs )
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
    sample_taxa_tsv = ch_taxonomy_tsv // channel: [ val(meta), tsv ]
    versions        = ch_versions     // channel: [ versions.yml ]
}

