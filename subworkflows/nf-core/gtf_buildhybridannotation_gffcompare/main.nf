include { GFFCOMPARE                      } from '../../../modules/nf-core/gffcompare'
include { GAWK as GAWK_FILTER             } from '../../../modules/nf-core/gawk'
include { BEDTOOLS_INTERSECT              } from '../../../modules/nf-core/bedtools/intersect'
include { GAWK as GAWK_CONCAT             } from '../../../modules/nf-core/gawk'

workflow GTF_BUILDHYBRIDANNOTATION_GFFCOMPARE {

    take:
    ch_novel_gtf         // channel: [ val(meta), path(novel_gtf) ]
    ch_reference_gtf     // channel: [ val(meta), path(reference_gtf) ]
    ch_backbone_gtf      // channel: [ val(meta), path(backbone_gtf) ]
    val_class_codes      // channel: val (String or List of gffcompare class codes to keep, e.g. "u" or "u,j,i")
    ch_blacklist_bed     // channel: [ val(meta), path(blacklist_bed) ] or Channel.empty() to skip the intersect

    main:

    // Classify every novel transcript against the reference annotation; this
    // tags transcripts in the annotated GTF with `class_code` attributes
    // (e.g. "u" intergenic, "j" alternative junction, "=" exact match).
    GFFCOMPARE(
        ch_novel_gtf,
        [[:], [], []],
        ch_reference_gtf
    )

    // Fold the requested class codes onto meta so GAWK_FILTER's closure-form
    // `ext.args` can read them at task launch via `meta.class_codes`.
    ch_filter_input = GFFCOMPARE.out.annotated_gtf
        .combine(val_class_codes)
        .map { meta, gtf, codes ->
            def codes_csv = codes instanceof List ? codes.join(',') : "${codes}"
            [ meta + [class_codes: codes_csv], gtf ]
        }

    // Drop transcripts whose class_code is not in the caller-supplied set
    // (also strips comment lines and rows with `.` strand).
    GAWK_FILTER(ch_filter_input, [], false)

    // Strip the transient `class_codes` key now that GAWK_FILTER has consumed
    // it, so downstream meta matches the caller's original shape and remains
    // joinable against other channels keyed on the same meta.
    ch_filtered = GAWK_FILTER.out.output
        .map { meta, gtf -> [ meta.findAll { k, _v -> k != 'class_codes' }, gtf ] }

    // Optional blacklist intersect: route surviving transcripts through
    // BEDTOOLS_INTERSECT only when the caller passes a blacklist BED. The
    // bedtools args (-v for negation, -s for strand, etc.) are caller policy
    // and live in the consumer pipeline's modules.config.
    ch_after_filter = ch_filtered
        .combine(ch_blacklist_bed.ifEmpty([[:], []]))
        .branch { _meta, _gtf, _bed_meta, bed ->
            with_blacklist: bed
            passthrough:    !bed
        }

    BEDTOOLS_INTERSECT(
        ch_after_filter.with_blacklist.map { meta, gtf, _bed_meta, bed -> [ meta, gtf, bed ] },
        [[:], []]
    )

    // Reunite the two branches into a single channel of survivors.
    ch_post_blacklist = BEDTOOLS_INTERSECT.out.intersect
        .mix(ch_after_filter.passthrough.map { meta, gtf, _bed_meta, _bed -> [ meta, gtf ] })

    // Pair the survivors with the backbone annotation and emit `[backbone, novel]`
    // (canonical first, so the concat awk sees backbone gene rows before any
    // novel transcript and can decide which gene_ids still need a synthetic row).
    ch_concat_input = ch_post_blacklist
        .combine(ch_backbone_gtf)
        .map { meta, novel_gtf, _backbone_meta, backbone_gtf -> [ meta, [ backbone_gtf, novel_gtf ] ] }

    // Concatenate backbone + novel, strip comments, and synthesise a `gene`
    // row for any novel gene_id absent from the backbone (using the first
    // transcript's coords as the gene span).
    GAWK_CONCAT(ch_concat_input, [], false)

    emit:
    hybrid_gtf          = GAWK_CONCAT.out.output       // channel: [ val(meta), path(hybrid.gtf) ]
    gffcompare_stats    = GFFCOMPARE.out.stats         // channel: [ val(meta), path(stats) ]
    gffcompare_tracking = GFFCOMPARE.out.tracking      // channel: [ val(meta), path(tracking) ]
    gffcompare_loci     = GFFCOMPARE.out.loci          // channel: [ val(meta), path(loci) ]
}
