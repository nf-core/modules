include { FREYJA_VARIANTS   } from '../../../modules/nf-core/freyja/variants'
include { FREYJA_UPDATE     } from '../../../modules/nf-core/freyja/update'
include { FREYJA_DEMIX      } from '../../../modules/nf-core/freyja/demix'
include { FREYJA_BOOT       } from '../../../modules/nf-core/freyja/boot'

workflow BAM_VARIANT_DEMIX_BOOT_FREYJA {

    take:
    ch_bam               // channel: [ val(meta), path(bam) ]
    ch_fasta             // channel: [ path(fasta) ]
    val_skip_boot        // value skip_boot
    val_repeats          // value repeats
    val_db_name          // string db_name
    ch_barcodes          // channel:  [ path(barcodes)]
    ch_lineages_meta     // channel:  [ path(lineages_meta)]
    ch_lineages_topology // channel:  [ path(lineages_topology)]

    main:
    //
    // Variant calling
    //
    FREYJA_VARIANTS (
        ch_bam,
        ch_fasta
    )
    ch_freyja_variants = FREYJA_VARIANTS.out.variants

    //
    // Update the database if none are given.
    //
    if (!ch_barcodes ) {
        // Throw warning if meta or topology are provided but no barcodes
        if (ch_lineages_meta || ch_lineages_topology) {
            log.warn("Lineage metadata or topology provided without barcodes. The database will be updated but the provided metadata/topology will be ignored.")
        }
        FREYJA_UPDATE (
            val_db_name
        )
        ch_barcodes          = FREYJA_UPDATE.out.barcodes
        ch_lineages_meta     = FREYJA_UPDATE.out.lineages_meta
        ch_lineages_topology = FREYJA_UPDATE.out.lineages_topology
    }

    //
    // demix and define minimum variant abundances
    //
    FREYJA_DEMIX (
        ch_freyja_variants,
        ch_barcodes,
        ch_lineages_meta,
        ch_lineages_topology
    )

    //
    // Perform bootstrapping to get more accurate estimates of abundances
    //
    ch_lineages   = channel.empty()
    ch_summarized = channel.empty()
    if (!val_skip_boot){
        FREYJA_BOOT (
            ch_freyja_variants,
            val_repeats,
            ch_barcodes,
            ch_lineages_meta,
            ch_lineages_topology
        )
        ch_lineages   = FREYJA_BOOT.out.lineages
        ch_summarized = FREYJA_BOOT.out.summarized
    }

    emit:
    variants          = FREYJA_VARIANTS.out.variants  // channel: [ val(meta), path(variants_tsv), path(depths_tsv) ]
    demix             = FREYJA_DEMIX.out.demix        // channel: [ val(meta), path(demix_tsv) ]
    lineages          = ch_lineages                   // channel: [ val(meta), path(lineages_csv) ]
    summarized        = ch_summarized                 // channel: [ val(meta), path(summarized_csv) ]
    barcodes          = ch_barcodes                   // channel: [ path(barcodes) ]
    lineages_meta     = ch_lineages_meta              // channel: [ path(lineages_meta) ]
    lineages_topology = ch_lineages_topology          // channel: [ path(lineages_topology) ]
    }
