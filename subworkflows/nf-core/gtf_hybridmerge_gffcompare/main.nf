include { GFFCOMPARE                      } from '../../../modules/nf-core/gffcompare'
include { GAWK as GAWK_FILTER             } from '../../../modules/nf-core/gawk'
include { BEDTOOLS_INTERSECT              } from '../../../modules/nf-core/bedtools/intersect'
include { GAWK as GAWK_CONCAT             } from '../../../modules/nf-core/gawk'

workflow GTF_HYBRIDMERGE_GFFCOMPARE {

    take:
    ch_novel_gtf         // channel: [ val(meta), path(novel_gtf) ]
    ch_reference_gtf     // channel: [ val(meta), path(reference_gtf) ]
    ch_backbone_gtf      // channel: [ val(meta), path(backbone_gtf) ]
    val_class_codes      // channel: val (String or List of gffcompare class codes to keep, e.g. "u" or "u,j,i")
    ch_blacklist_bed     // channel: [ val(meta), path(blacklist_bed) ] or channel.empty() to skip the intersect

    main:

    // Bundled awk programs (see ./awk/{filter,concat}.awk); shipped alongside
    // the subworkflow and passed to GAWK as `program_file` so the awk source
    // stays readable as awk rather than as an escaped Groovy string.
    ch_filter_awk = file("${moduleDir}/awk/filter.awk", checkIfExists: true)
    ch_concat_awk = file("${moduleDir}/awk/concat.awk", checkIfExists: true)

    // Tag every novel transcript with a gffcompare class_code attribute.
    GFFCOMPARE(
        ch_novel_gtf,
        [[:], [], []],
        ch_reference_gtf
    )

    // Fold class codes onto meta so GAWK_FILTER's ext.args closure can read them at task launch.
    ch_filter_input = GFFCOMPARE.out.annotated_gtf
        .combine(val_class_codes)
        .map { meta, gtf, codes ->
            def codes_csv = codes instanceof List ? codes.join(',') : "${codes}"
            [ meta + [class_codes: codes_csv], gtf ]
        }

    // Drop transcripts whose class_code is not in the caller-supplied set.
    GAWK_FILTER(ch_filter_input, ch_filter_awk, false)

    // Branch: send survivors through BEDTOOLS_INTERSECT only when a blacklist BED was supplied.
    // Bedtools args (-v for negation, -s for strand) are caller policy in modules.config.
    ch_after_filter = GAWK_FILTER.out.output
        .combine(ch_blacklist_bed.ifEmpty([[:], []]))
        .branch { _meta, _gtf, _bed_meta, bed ->
            with_blacklist: bed
            passthrough:    !bed
        }

    BEDTOOLS_INTERSECT(
        ch_after_filter.with_blacklist.map { meta, gtf, _bed_meta, bed -> [ meta, gtf, bed ] },
        [[:], []]
    )

    ch_post_blacklist = BEDTOOLS_INTERSECT.out.intersect
        .mix(ch_after_filter.passthrough.map { meta, gtf, _bed_meta, _bed -> [ meta, gtf ] })

    // Adopt the backbone's meta (the hybrid identity follows the backbone) and stage
    // [backbone, novel] in that order so backbone gene rows reach the concat awk first.
    ch_concat_input = ch_post_blacklist
        .combine(ch_backbone_gtf)
        .map { _novel_meta, novel_gtf, backbone_meta, backbone_gtf ->
            [ backbone_meta, [ backbone_gtf, novel_gtf ] ]
        }

    // Concatenate and synthesise missing gene rows with union spans.
    GAWK_CONCAT(ch_concat_input, ch_concat_awk, false)

    emit:
    hybrid_gtf          = GAWK_CONCAT.out.output       // channel: [ val(meta), path(hybrid.gtf) ]
    gffcompare_stats    = GFFCOMPARE.out.stats         // channel: [ val(meta), path(stats) ]
    gffcompare_tracking = GFFCOMPARE.out.tracking      // channel: [ val(meta), path(tracking) ]
    gffcompare_loci     = GFFCOMPARE.out.loci          // channel: [ val(meta), path(loci) ]
}
