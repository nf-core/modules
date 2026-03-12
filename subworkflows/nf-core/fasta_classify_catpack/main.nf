/*
 * CAT/BAT/RAT: tools for taxonomic classification of contigs and metagenome-assembled genomes (MAGs)
 */
include { CATPACK_ADDNAMES as CATPACK_ADDNAMES_BINS     } from '../../../modules/nf-core/catpack/addnames/main'
include { CATPACK_ADDNAMES as CATPACK_ADDNAMES_UNBINS   } from '../../../modules/nf-core/catpack/addnames/main'
include { CATPACK_BINS                                  } from '../../../modules/nf-core/catpack/bins/main'
include { CATPACK_CONTIGS as CATPACK_UNBINS             } from '../../../modules/nf-core/catpack/contigs/main'
include { CATPACK_DOWNLOAD                              } from '../../../modules/nf-core/catpack/download/main'
include { CATPACK_PREPARE                               } from '../../../modules/nf-core/catpack/prepare/main'
include { CATPACK_SUMMARISE as CATPACK_SUMMARISE_BINS   } from '../../../modules/nf-core/catpack/summarise/main'
include { CATPACK_SUMMARISE as CATPACK_SUMMARISE_UNBINS } from '../../../modules/nf-core/catpack/summarise/main'
include { UNTAR as CAT_DB_UNTAR                         } from '../../../modules/nf-core/untar/main'

workflow FASTA_CLASSIFY_CATPACK {

    take:
    ch_bins                   // channel: [ val(meta), path(fasta) ] - binned MAGs/contigs
    ch_unbins                 // channel: [ val(meta), path(fasta) ] - unbinned contigs
    ch_cat_db                 // channel: [ val(meta), path(db) ] - pre-built db as directory (with db/ and tax/ subdirs) or .tar.gz
                              //          provide channel.empty() to trigger automatic download via ch_cat_db_download_id
    ch_cat_db_download_id     // channel: [ val(meta), val(db_id) ] - db ID for CATPACK_DOWNLOAD (e.g. 'nr')
                              //          provide channel.empty() if supplying a pre-built db via ch_cat_db
    classify_unbinned         // val: boolean - whether to classify unbinned contigs
    allow_unofficial_lineages // val: boolean - whether to allow unofficial lineages in the summarise step
    bin_suffix                // val: string - file extension of bin FASTA files (e.g. '.fa')

    main:
    ch_versions = channel.empty()

    //
    // Database preparation
    //

    // Handle pre-built db: untar if compressed, or use directory directly
    ch_cat_db
        .branch { _meta, db ->
            tar: db.name.endsWith('.tar.gz')
            dir: db.isDirectory()
            other: true
         }
        .set { ch_cat_db_input }

ch_cat_db.other.subscribe { db -> 
    exit("Error: A DB was provided to FASTA_CLASSIFY_CATPACK that is not a `.tar.gz` or a directory!")
}
    CAT_DB_UNTAR(ch_cat_db_input.tar)
    ch_versions = ch_versions.mix(CAT_DB_UNTAR.out.versions_untar.first())

    ch_prepared_from_dir = ch_cat_db_input.dir
        .mix(CAT_DB_UNTAR.out.untar)
        .multiMap { meta, dir ->
            db:       [meta, dir / 'db']
            taxonomy: [meta, dir / 'tax']
        }

    // Download and prepare db from scratch if no pre-built db provided
    CATPACK_DOWNLOAD(ch_cat_db_download_id)
    ch_versions = ch_versions.mix(CATPACK_DOWNLOAD.out.versions_catpack.first())

    CATPACK_PREPARE(
        CATPACK_DOWNLOAD.out.fasta,
        CATPACK_DOWNLOAD.out.names.map   { _meta, names -> names },
        CATPACK_DOWNLOAD.out.nodes.map   { _meta, nodes -> nodes },
        CATPACK_DOWNLOAD.out.acc2tax.map { _meta, acc2tax -> acc2tax },
    )
    ch_versions = ch_versions.mix(CATPACK_PREPARE.out.versions_catpack.first())

    // Combine db sources - one of these channels will be empty depending on inputs
    ch_db       = ch_prepared_from_dir.db.mix(CATPACK_PREPARE.out.db).first()
    ch_taxonomy = ch_prepared_from_dir.taxonomy.mix(CATPACK_PREPARE.out.taxonomy).first()

    //
    // Bin taxonomic classification
    //

    CATPACK_BINS(
        ch_bins,
        ch_db,
        ch_taxonomy,
        [[:], []],
        [[:], []],
        bin_suffix,
    )
    ch_versions = ch_versions.mix(CATPACK_BINS.out.versions_catpack.first())

    CATPACK_ADDNAMES_BINS(CATPACK_BINS.out.bin2classification, ch_taxonomy)
    ch_versions = ch_versions.mix(CATPACK_ADDNAMES_BINS.out.versions_catpack.first())

    ch_summarise_bins = channel.empty()
    if (!allow_unofficial_lineages) {
        CATPACK_SUMMARISE_BINS(CATPACK_ADDNAMES_BINS.out.txt, [[:], []])
        ch_versions = ch_versions.mix(CATPACK_SUMMARISE_BINS.out.versions_catpack.first())
        ch_summarise_bins = CATPACK_SUMMARISE_BINS.out.txt
    }

    //
    // Unbinned contigs taxonomic classification (optional)
    //

    ch_unbins_addnames  = channel.empty()
    ch_summarise_unbins = channel.empty()

    if (classify_unbinned) {
        CATPACK_UNBINS(
            ch_unbins,
            ch_db,
            ch_taxonomy,
            [[:], []],
            [[:], []],
        )
        ch_versions = ch_versions.mix(CATPACK_UNBINS.out.versions_catpack.first())

        CATPACK_ADDNAMES_UNBINS(CATPACK_UNBINS.out.contig2classification, ch_taxonomy)
        ch_versions = ch_versions.mix(CATPACK_ADDNAMES_UNBINS.out.versions_catpack.first())

        ch_unbins_addnames = CATPACK_ADDNAMES_UNBINS.out.txt

        if (!allow_unofficial_lineages) {
            ch_unbin_input = CATPACK_ADDNAMES_UNBINS.out.txt
                .join(ch_unbins)
                .multiMap { meta, names, contigs ->
                    names:   [meta, names]
                    contigs: [meta, contigs]
                }

            CATPACK_SUMMARISE_UNBINS(ch_unbin_input.names, ch_unbin_input.contigs)
            ch_versions = ch_versions.mix(CATPACK_SUMMARISE_UNBINS.out.versions_catpack.first())
            ch_summarise_unbins = CATPACK_SUMMARISE_UNBINS.out.txt
        }
    }

    emit:
    bat_classification      = CATPACK_ADDNAMES_BINS.out.txt // channel: [ val(meta), path(txt) ]
    bat_summary             = ch_summarise_bins             // channel: [ val(meta), path(txt) ]
    unbinned_classification = ch_unbins_addnames            // channel: [ val(meta), path(txt) ]
    unbinned_summary        = ch_summarise_unbins           // channel: [ val(meta), path(txt) ]
    versions                = ch_versions                   // channel: versions
}
