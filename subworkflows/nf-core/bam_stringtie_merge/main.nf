include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'


workflow BAM_STRINGTIE_MERGE {
    take:
    bam_sorted // channel: [ meta, bam ]
    ch_chrgtf // channel: [ meta, gtf ]

    main:
    ch_stringtie_gtfs = channel.empty()

    STRINGTIE_STRINGTIE(
        bam_sorted,
        ch_chrgtf.map { _meta, gtf -> [gtf] },
    )

    STRINGTIE_STRINGTIE.out.transcript_gtf.set { stringtie_gtfs }

    STRINGTIE_MERGE(
        stringtie_gtfs,
        ch_chrgtf
    )
    ch_stringtie_gtfs = STRINGTIE_MERGE.out.merged_gtf

    emit:
    stringtie_gtf = ch_stringtie_gtfs // channel: [ meta, gtf ]
}
