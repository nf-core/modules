//
// Short read sequencing data QC using different tools
//
include { FASTQC       } from '../../../modules/nf-core/fastqc/main'
include { SEQFU_CHECK  } from '../../../modules/nf-core/seqfu/check/main'
include { SEQFU_STATS  } from '../../../modules/nf-core/seqfu/stats/main'
include { SEQKIT_STATS } from '../../../modules/nf-core/seqkit/stats/main'
include { SEQTK_COMP   } from '../../../modules/nf-core/seqtk/comp/main'

workflow FASTQ_GENERATE_STATISTICS {

    take:
    ch_reads          // channel: [ val(meta), [ fastq ] ]
    skip_fastqc       // boolean
    skip_seqfu_check  // boolean
    skip_seqfu_stats  // boolean
    skip_seqkit_stats // boolean
    skip_seqtk_comp   // boolean

    main:

    ch_versions = Channel.empty()

    if (!skip_fastqc) {
        FASTQC( ch_reads )
    }

    if (!skip_seqfu_check){
        SEQFU_CHECK( ch_reads )
        ch_versions = ch_versions.mix(SEQFU_CHECK.out.versions.first())
    }

    if (!skip_seqfu_stats) {
        SEQFU_STATS ( ch_reads )
        ch_versions = ch_versions.mix(SEQFU_STATS.out.versions.first())
    }

    if (!skip_seqkit_stats) {
        SEQKIT_STATS ( ch_reads )
        ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions.first())
    }

    if (!skip_seqtk_comp) {
        SEQTK_COMP ( ch_reads )
        ch_versions = ch_versions.mix(SEQTK_COMP.out.versions.first())
    }

    emit:
    fastqc_html   = FASTQC.out.html
    fastqc_zip    = FASTQC.out.zip
    seqfu_check   = SEQFU_CHECK.out.check
    seqfu_stats   = SEQFU_STATS.out.stats
    seqfu_multiqc = SEQFU_STATS.out.multiqc
    seqkit_stats  = SEQKIT_STATS.out.stats
    seqtk_stats   = SEQTK_COMP.out.seqtk_stats
    versions      = ch_versions
}
