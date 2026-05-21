include { FASTQC       } from '../../../modules/nf-core/fastqc'
include { SEQFU_CHECK  } from '../../../modules/nf-core/seqfu/check'
include { SEQFU_STATS  } from '../../../modules/nf-core/seqfu/stats'
include { SEQKIT_STATS } from '../../../modules/nf-core/seqkit/stats'
include { SEQTK_COMP   } from '../../../modules/nf-core/seqtk/comp'

workflow FASTQ_QC_STATS {
    take:
    ch_reads // channel: [ val(meta), [ fastq ] ]
    skip_fastqc // boolean
    skip_seqfu_check // boolean
    skip_seqfu_stats // boolean
    skip_seqkit_stats // boolean
    skip_seqtk_comp // boolean

    main:
    ch_fastqc_html = channel.empty()
    ch_fastqc_zip = channel.empty()
    ch_seqfu_check = channel.empty()
    ch_seqfu_stats = channel.empty()
    ch_seqfu_multiqc = channel.empty()
    ch_seqkit_stats = channel.empty()
    ch_seqtk_stats = channel.empty()

    if (!skip_fastqc) {
        FASTQC(ch_reads)
        ch_fastqc_html = FASTQC.out.html
        ch_fastqc_zip = FASTQC.out.zip
    }

    if (!skip_seqfu_check) {
        SEQFU_CHECK(ch_reads)
        ch_seqfu_check = SEQFU_CHECK.out.check
    }

    if (!skip_seqfu_stats) {
        SEQFU_STATS(ch_reads)
        ch_seqfu_stats = SEQFU_STATS.out.stats
        ch_seqfu_multiqc = SEQFU_STATS.out.multiqc
    }

    if (!skip_seqkit_stats) {
        SEQKIT_STATS(ch_reads)
        ch_seqkit_stats = SEQKIT_STATS.out.stats
    }

    if (!skip_seqtk_comp) {
        SEQTK_COMP(ch_reads)
        ch_seqtk_stats = SEQTK_COMP.out.seqtk_stats
    }

    emit:
    fastqc_html   = ch_fastqc_html
    fastqc_zip    = ch_fastqc_zip
    seqfu_check   = ch_seqfu_check
    seqfu_stats   = ch_seqfu_stats
    seqfu_multiqc = ch_seqfu_multiqc
    seqkit_stats  = ch_seqkit_stats
    seqtk_stats   = ch_seqtk_stats
}
