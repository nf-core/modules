include { NGSCHECKMATE_FASTQ } from '../../../modules/nf-core/ngscheckmate/fastq/main'
include { NGSCHECKMATE_VAFNCM } from '../../../modules/nf-core/ngscheckmate/vafncm/main'

workflow FASTQ_NGSCHECKMATE {

    take:
    ch_fastq            // channel: [ val(meta1), fastq ]
    ch_snp_pt           // channel: [ val(meta2), snp_pt ]

    main:

    ch_versions = Channel.empty()

    NGSCHECKMATE_FASTQ (ch_fastq, ch_snp_pt )
    ch_versions = ch_versions.mix(NGSCHECKMATE_FASTQ.out.versions.first())

    NGSCHECKMATE_FASTQ
        .out
        .vaf
        .map{meta, vaf -> vaf}    // discard individual metas
        .collect()                // group into one channel
        .map{files -> [files]}    // make the channel into [vaf1, vaf2, ...]
        .set {ch_collected_vafs}

    ch_snp_pt
        .map{meta, snp_pt -> meta} // use the snp_bed file meta as the meta for the merged channel
        .combine(ch_collected_vafs) // add the vaf files after the meta, now looks like [meta, [vaf1, vaf2, ... ] ]
        .set {ch_vafs}

    NGSCHECKMATE_VAFNCM (ch_vafs)
    ch_versions = ch_versions.mix(NGSCHECKMATE_VAFNCM.out.versions.first())

    emit:

    corr_matrix  = NGSCHECKMATE_VAFNCM.out.corr_matrix  // channel: [ meta, corr_matrix ]
    matched      = NGSCHECKMATE_VAFNCM.out.matched      // channel: [ meta, matched ]
    all          = NGSCHECKMATE_VAFNCM.out.all          // channel: [ meta, all ]
    vaf          = NGSCHECKMATE_FASTQ.out.vaf           // channel: [ meta, vaf ]
    pdf          = NGSCHECKMATE_VAFNCM.out.pdf          // channel: [ meta, pdf ], optional
    versions     = ch_versions                          // channel: [ versions.yml ]
}

