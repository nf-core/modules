include { FREYJA_VARIANTS   } from '../../../modules/nf-core/freyja/variants'
include { FREYJA_UPDATE     } from '../../../modules/nf-core/freyja/update'
include { FREYJA_DEMIX      } from '../../../modules/nf-core/freyja/demix'
include { FREYJA_BOOT       } from '../../../modules/nf-core/freyja/boot'

workflow BAM_VARIANT_DEMIX_BOOT_FREYJA {

    take:
    ch_bam          // channel: [ val(meta), [ /path/to/bam ] ]
    fasta           // channel: [ val(meta), [path/to/reference.fasta] ]
    repeats         // value repeats
    db_name         // string db_name
    barcodes        // /path/to/barcodes.csv
    lineage_meta    // /path/to/lineage_metadata

    main:

    ch_versions = Channel.empty()

    //
    // Variant calling
    //
    ch_freyja_variants = Channel.empty()
    ch_freyja_depths    = Channel.empty()
    FREYJA_VARIANTS(
        ch_bam,
        fasta)
    ch_freyja_variants = FREYJA_VARIANTS.out.variants
    ch_freyja_depths   = FREYJA_VARIANTS.out.depths

    ch_versions=ch_versions.mix(FREYJA_VARIANTS.out.versions)

    //
    // Update the database if none are given.
    //
    if(barcodes==[] || lineage_meta == []){
        FREYJA_UPDATE(db_name)
        barcodes=FREYJA_UPDATE.out.barcodes
        lineage_meta=FREYJA_UPDATE.out.lineages_meta
        ch_versions=ch_versions.mix(FREYJA_UPDATE.out.versions)
    }

    //
    // demix and define minimum variant abundances
    //
    ch_freyja_abundacies    = Channel.empty()
    FREYJA_DEMIX(
        ch_freyja_variants,
        ch_freyja_depths,
        barcodes,
        lineage_meta)
    ch_freyja_demix = FREYJA_DEMIX.out.demix
    ch_versions=ch_versions.mix(FREYJA_DEMIX.out.versions)

    //
    // Perform bootstrapping to get more accurate estimates of abundancies
    //
    FREYJA_BOOT(
        ch_freyja_variants,
        ch_freyja_depths,
        repeats,
        barcodes,
        lineage_meta)

    ch_versions=ch_versions.mix(FREYJA_BOOT.out.versions)

    emit:
    variants    = FREYJA_VARIANTS.out.variants  // channel: [ val(meta), [ variants.tsv ] ]
    depths      = FREYJA_VARIANTS.out.depths    // channel: [ val(meta), [ depths.tsv ] ]
    demix       = FREYJA_DEMIX.out.demix        // channel: [ val(meta), [ demix.tsv ] ]
    lineages    = FREYJA_BOOT.out.lineages      // channel: [ val(meta), [ lineages.csv ] ]
    summarized  = FREYJA_BOOT.out.summarized    // channel: [ val(meta), [ summarized.csv ] ]
    versions    = ch_versions                   // channel: [ versions.yml ]
}

