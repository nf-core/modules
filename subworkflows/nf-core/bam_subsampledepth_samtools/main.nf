include { SAMTOOLS_DEPTH } from '../../../modules/nf-core/samtools/depth'
include { GAWK           } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view'

workflow BAM_SUBSAMPLEDEPTH_SAMTOOLS {

    take:
    ch_bam_bai_depth    // channel: [ val(meta), path(bam), path(bai), val(depth) ]
    ch_fasta            // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions      = Channel.empty()

    // Compute mean depth
    SAMTOOLS_DEPTH(ch_bam_bai_depth.map{ it[0..2] }, [[], []])
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    // Use GAWK to get mean depth
    GAWK(SAMTOOLS_DEPTH.out.tsv, [])
    ch_versions = ch_versions.mix(GAWK.out.versions.first())

    // Compute downsampling factor
    ch_mean_depth = GAWK.out.output
        .splitCsv(header: false, sep:'\t')
        .map{ meta, row ->
            [ meta, row[0] as Float ]
        }

    // Add all necessary channel for downsampling
    ch_input_subsample = ch_bam_bai_depth
        .join(ch_mean_depth)
        .map{ meta, bam, index, depth, mean ->
            [ meta + ['subsample_fraction': depth as Float / mean, 'depth': depth ], bam, index ]
        }

    // Downsample
    SAMTOOLS_VIEW(
        ch_input_subsample,
        ch_fasta,
        [],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Aggregate bam and index
    ch_bam_subsampled = SAMTOOLS_VIEW.out.bam.mix(SAMTOOLS_VIEW.out.cram, SAMTOOLS_VIEW.out.sam)
        .join(SAMTOOLS_VIEW.out.bai.mix(SAMTOOLS_VIEW.out.crai, SAMTOOLS_VIEW.out.csi))

    emit:
    bam_subsampled    = ch_bam_subsampled             // channel: [ val(meta), path(bam), path(csi) ]
    versions          = ch_versions                   // channel: [ path(versions.yml) ]
}
