include { WISECONDORX_CONVERT } from '../../../modules/nf-core/wisecondorx/convert/main'
include { WISECONDORX_PREDICT } from '../../../modules/nf-core/wisecondorx/predict/main'

workflow BAM_CNV_WISECONDORX {

    take:
    ch_bam          // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta        // channel: [ val(meta2), path(fasta) ]
    ch_fai          // channel: [ val(meta3), path(fai) ]
    ch_ref          // channel: [ val(meta4), path(reference) ]
    ch_blacklist    // channel: [ val(meta5), path(blacklist) ]

    main:

    ch_versions = Channel.empty()

    WISECONDORX_CONVERT(
        ch_bam,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(WISECONDORX_CONVERT.out.versions.first())

    WISECONDORX_PREDICT(
        WISECONDORX_CONVERT.out.npz,
        ch_ref,
        ch_blacklist
    )
    ch_versions = ch_versions.mix(WISECONDORX_PREDICT.out.versions.first())

    emit:
    aberrations_bed = WISECONDORX_PREDICT.out.aberrations_bed   // channel: [ val(meta), path(bed) ]
    bins_bed        = WISECONDORX_PREDICT.out.bins_bed          // channel: [ val(meta), path(bed) ]
    segments_bed    = WISECONDORX_PREDICT.out.segments_bed      // channel: [ val(meta), path(bed) ]
    chr_statistics  = WISECONDORX_PREDICT.out.chr_statistics    // channel: [ val(meta), path(txt) ]
    chr_plots       = WISECONDORX_PREDICT.out.chr_plots         // channel: [ val(meta), [ path(png), path(png), ... ] ]
    genome_plot     = WISECONDORX_PREDICT.out.genome_plot       // channel: [ val(meta), path(png) ]

    versions        = ch_versions                               // channel: path(versions.yml)
}
