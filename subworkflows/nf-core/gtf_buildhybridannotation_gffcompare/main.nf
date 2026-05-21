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

    GFFCOMPARE(
        ch_novel_gtf,
        [[:], [], []],
        ch_reference_gtf
    )

    ch_filter_input = GFFCOMPARE.out.annotated_gtf
        .combine(val_class_codes)
        .map { meta, gtf, codes ->
            def codes_csv = codes instanceof List ? codes.join(',') : "${codes}"
            [ meta + [class_codes: codes_csv], gtf ]
        }

    GAWK_FILTER(ch_filter_input, [], false)

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

    ch_concat_input = ch_post_blacklist
        .combine(ch_backbone_gtf)
        .map { meta, novel_gtf, _backbone_meta, backbone_gtf -> [ meta, [ backbone_gtf, novel_gtf ] ] }

    GAWK_CONCAT(ch_concat_input, [], false)

    emit:
    hybrid_gtf          = GAWK_CONCAT.out.output       // channel: [ val(meta), path(hybrid.gtf) ]
    gffcompare_stats    = GFFCOMPARE.out.stats         // channel: [ val(meta), path(stats) ]
    gffcompare_tracking = GFFCOMPARE.out.tracking      // channel: [ val(meta), path(tracking) ]
    gffcompare_loci     = GFFCOMPARE.out.loci          // channel: [ val(meta), path(loci) ]
}
