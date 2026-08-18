/*
 * CAT/BAT/RAT: tools for taxonomic classification of contigs and metagenome-assembled genomes (MAGs)
 */
include { CATPACK_ADDNAMES as CATPACK_ADDNAMES_BINS    } from '../../../modules/nf-core/catpack/addnames/main'
include { CATPACK_ADDNAMES as CATPACK_ADDNAMES_CONTIGS } from '../../../modules/nf-core/catpack/addnames/main'
include { CATPACK_BINS                                 } from '../../../modules/nf-core/catpack/bins/main'
include { CATPACK_CONTIGS                              } from '../../../modules/nf-core/catpack/contigs/main'
include { CATPACK_DOWNLOAD                             } from '../../../modules/nf-core/catpack/download/main'
include { CATPACK_PREPARE                              } from '../../../modules/nf-core/catpack/prepare/main'
include { CATPACK_SUMMARISE as CATPACK_SUMMARISE_BINS    } from '../../../modules/nf-core/catpack/summarise/main'
include { CATPACK_SUMMARISE as CATPACK_SUMMARISE_CONTIGS } from '../../../modules/nf-core/catpack/summarise/main'
include { UNTAR as CAT_DB_UNTAR                        } from '../../../modules/nf-core/untar/main'

workflow FASTA_CLASSIFY_CATPACK {

    take:
    ch_bins               // channel: [ val(meta), path(fasta) ] - binned MAGs/contigs
    ch_contigs            // channel: [ val(meta), path(fasta) ] - contigs; provide channel.empty() to skip contig classification
    ch_cat_db             // channel: [ val(meta), path(db) ] - pre-built db as directory (with db/ and tax/ subdirs) or .tar.gz
                          //          provide channel.empty() to trigger automatic download via ch_cat_db_download_id
    ch_cat_db_download_id // channel: [ val(meta), val(db_id) ] - db ID for CATPACK_DOWNLOAD (e.g. 'nr')
                          //          provide channel.empty() if supplying a pre-built db via ch_cat_db
                          //          supplying both ch_cat_db and ch_cat_db_download_id will cause a runtime error
    run_summarise         // val: boolean - whether to run CATPACK_SUMMARISE; requires ext.args = "--only_official" on CATPACK_ADDNAMES_BINS/CONTIGS
    bin_suffix            // val: string - file extension of bin FASTA files (e.g. '.fa' or '.fasta')

    main:

    //
    // Database preparation
    //

    // Handle pre-built db: untar if compressed, or use directory directly
    ch_cat_db_input = ch_cat_db
        .branch { _meta, db ->
            tar:   db.name.endsWith('.tar.gz')
            dir:   db.isDirectory()
            other: true
         }

    ch_cat_db_input.other.subscribe { _meta, _db ->
        error("Error: A DB was provided to FASTA_CLASSIFY_CATPACK that is not a `.tar.gz` or a directory.")
    }

    CAT_DB_UNTAR(ch_cat_db_input.tar)

    ch_prepared_from_dir = ch_cat_db_input.dir
        .mix(CAT_DB_UNTAR.out.untar)
        .multiMap { meta, dir ->
            db:       [meta, dir / 'db']
            taxonomy: [meta, dir / 'tax']
        }

    // Download and prepare db from scratch if no pre-built db provided
    CATPACK_DOWNLOAD(ch_cat_db_download_id)

    CATPACK_PREPARE(
        CATPACK_DOWNLOAD.out.fasta,
        CATPACK_DOWNLOAD.out.names.map   { _meta, names -> names },
        CATPACK_DOWNLOAD.out.nodes.map   { _meta, nodes -> nodes },
        CATPACK_DOWNLOAD.out.acc2tax.map { _meta, acc2tax -> acc2tax },
    )

    // Combine db sources - one of these channels will be empty depending on inputs
    // Guard: fail if both ch_cat_db and ch_cat_db_download_id are provided simultaneously.
    // .combine() only emits when both channels have at least one element.
    ch_prepared_from_dir.db.combine(CATPACK_PREPARE.out.db).subscribe {
        error("Error: Both a pre-built DB and a download ID were provided to FASTA_CLASSIFY_CATPACK! Provide only one via ch_cat_db or ch_cat_db_download_id.")
    }

    ch_db       = ch_prepared_from_dir.db.mix(CATPACK_PREPARE.out.db).first()
    ch_taxonomy = ch_prepared_from_dir.taxonomy.mix(CATPACK_PREPARE.out.taxonomy).first()

    //
    // Bin taxonomic classification (optional - skipped when ch_bins is channel.empty())
    //

    CATPACK_BINS(
        ch_bins,
        ch_db,
        ch_taxonomy,
        [[:], []],
        [[:], []],
        bin_suffix,
    )

    CATPACK_ADDNAMES_BINS(CATPACK_BINS.out.bin2classification, ch_taxonomy)

    ch_bat_summary = channel.empty()
    if (run_summarise) {
        CATPACK_SUMMARISE_BINS(CATPACK_ADDNAMES_BINS.out.txt, [[:], []])
        ch_bat_summary = CATPACK_SUMMARISE_BINS.out.txt
    }

    //
    // Contig taxonomic classification (optional - skipped when ch_contigs is channel.empty())
    //

    CATPACK_CONTIGS(
        ch_contigs,
        ch_db,
        ch_taxonomy,
        [[:], []],
        [[:], []],
    )

    CATPACK_ADDNAMES_CONTIGS(CATPACK_CONTIGS.out.contig2classification, ch_taxonomy)

    ch_contigs_summary = channel.empty()
    if (run_summarise) {
        ch_contigs_input = CATPACK_ADDNAMES_CONTIGS.out.txt
            .join(ch_contigs)
            .multiMap { meta, names, contigs ->
                names:   [meta, names]
                contigs: [meta, contigs]
            }

        CATPACK_SUMMARISE_CONTIGS(ch_contigs_input.names, ch_contigs_input.contigs)
        ch_contigs_summary = CATPACK_SUMMARISE_CONTIGS.out.txt
    }

    emit:
    bin2classification      = CATPACK_BINS.out.bin2classification        // channel: [ val(meta), path(txt) ]
    bat_classification      = CATPACK_ADDNAMES_BINS.out.txt              // channel: [ val(meta), path(txt) ]
    bat_summary             = ch_bat_summary                             // channel: [ val(meta), path(txt) ]
    contig2classification   = CATPACK_CONTIGS.out.contig2classification  // channel: [ val(meta), path(txt) ]
    contigs_classification  = CATPACK_ADDNAMES_CONTIGS.out.txt           // channel: [ val(meta), path(txt) ]
    contigs_summary         = ch_contigs_summary                         // channel: [ val(meta), path(txt) ]
}
