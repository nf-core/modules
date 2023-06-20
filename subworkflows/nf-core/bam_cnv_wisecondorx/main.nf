include { WISECONDORX_CONVERT } from '../../../modules/nf-core/wisecondorx/convert/main'
include { WISECONDORX_PREDICT } from '../../../modules/nf-core/wisecondorx/predict/main'

workflow  {

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
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

