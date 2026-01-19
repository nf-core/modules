include { NGSCHECKMATE_FASTQ  } from '../../../modules/nf-core/ngscheckmate/fastq/main'
include { NGSCHECKMATE_VAFNCM } from '../../../modules/nf-core/ngscheckmate/vafncm/main'

workflow FASTQ_NGSCHECKMATE {
    take:
    ch_fastq  // channel: [ val(meta1), fastq ]
    ch_snp_pt // channel: [ val(meta2), snp_pt ]

    main:
    NGSCHECKMATE_FASTQ(ch_fastq, ch_snp_pt.first())

    NGSCHECKMATE_FASTQ.out.vaf.map { _meta, vaf -> vaf }.collect().map { files -> [files] }.set { ch_collected_vafs }

    ch_snp_pt
        .map { meta, _snp_pt -> meta }
        .combine(ch_collected_vafs)
        .set { ch_vafs }

    NGSCHECKMATE_VAFNCM(ch_vafs)

    emit:
    corr_matrix = NGSCHECKMATE_VAFNCM.out.corr_matrix // channel: [ meta, corr_matrix ]
    matched     = NGSCHECKMATE_VAFNCM.out.matched     // channel: [ meta, matched ]
    all         = NGSCHECKMATE_VAFNCM.out.all         // channel: [ meta, all ]
    vaf         = NGSCHECKMATE_FASTQ.out.vaf          // channel: [ meta, vaf ]
    pdf         = NGSCHECKMATE_VAFNCM.out.pdf         // channel: [ meta, pdf ], optional
}
