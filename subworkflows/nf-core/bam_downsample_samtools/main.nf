include { SAMTOOLS_DEPTH } from '../../../modules/nf-core/samtools/depth'
include { GAWK           } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view'

workflow BAM_DOWNSAMPLE_SAMTOOLS {

    take:
    ch_bam    // channel: [ [id], bam, bai ]
    ch_fasta  // channel: [ [genome], fasta, fai ]
    depth     // channel: [ depth ]

    main:
    ch_versions      = Channel.empty()

    // Compute mean depth
    SAMTOOLS_DEPTH(ch_bam, [[], []])
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    // Use GAWK to get mean depth
    GAWK(SAMTOOLS_DEPTH.out.tsv, [])
    ch_versions = ch_versions.mix(GAWK.out.versions.first())

    // Compute downsampling factor
    ch_depth_factor = GAWK.out.output
        .splitCsv(header: false, sep:'\t')
        .map{ metaI, row ->
            [ metaI, row[0] as Float ]
        }
        .combine(depth)
        .map{ metaI, mean, depth ->
            [ metaI, metaI + ["depth": depth], depth as Float / mean ]
        }

    // Add all necessary channel for downsampling
    ch_input_downsample = ch_bam
        .combine(ch_depth_factor, by : 0)
        .map{ metaI, bam, index, metaID, depth ->
            [ metaID, bam, index, [], depth ]
        }

    // Downsample
    SAMTOOLS_VIEW(
        ch_input_downsample,
        ch_fasta
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Aggregate bam and index
    ch_bam_emul = SAMTOOLS_VIEW.out.bam
        .join(SAMTOOLS_VIEW.out.csi)

    emit:
    bam_downsampled   = ch_bam_emul                    // channel: [ [id, depth], bam, bai ]
    versions          = ch_versions                    // channel: [ versions.yml ]
}
