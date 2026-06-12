include { SAMTOOLS_COVERAGE } from '../../../modules/nf-core/samtools/coverage'
include { SAMTOOLS_VIEW     } from '../../../modules/nf-core/samtools/view'

workflow BAM_SUBSAMPLEDEPTH_SAMTOOLS {

    take:
    ch_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
    ch_depth      // channel: [ val(meta), val(depth)]
    ch_fasta      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_regions    // channel (optional): [ val(meta), path(region) ]

    main:

    ch_regions_branched = ch_regions
        .branch { _meta, region ->
            with_bed: region
            without_bed: true
        }

    ch_regions_parsed = ch_regions_branched.with_bed
        .splitCsv(header: false, sep:'\t')
        .map{ meta, row -> [
            meta, "${row[0]}:${row[1]}-${row[2]}"
        ]}
    // For samples WITH BED file: compute coverage per region
    ch_coverage_with_regions = ch_bam_bai
        .combine(ch_regions_parsed)
        .map{ meta, bam, bai, _metaR, region ->
            [ meta + ["region": region], bam, bai ]
        }

    // For samples WITHOUT BED file: compute whole-genome coverage
    ch_coverage_without_regions = ch_bam_bai
        .combine(ch_regions_branched.without_bed)
        .map{ meta, bam, bai, _metaR, _region_file -> [ meta, bam, bai ] }

    SAMTOOLS_COVERAGE(
        ch_coverage_with_regions
            .mix(ch_coverage_without_regions),
        ch_fasta
    )

    // Compute downsampling factor
    ch_mean_depth = SAMTOOLS_COVERAGE.out.coverage
        .splitCsv(header: true, sep:'\t')
        .filter{ _meta, row -> row.meandepth as Float > 0 }
        .map{ meta, row ->
            def keys_to_keep = meta.keySet() - ['region']
            [
                meta.subMap(keys_to_keep),
                row.meandepth as Float
            ]
        }
        .groupTuple()
        .map{ meta, mean ->
            [ meta, mean.sum() / mean.size() ]
        }

    ch_input_subsample = ch_bam_bai
        .join(ch_mean_depth)
        .combine(ch_depth)
        .map{ meta, bam, index, mean, metaD, depth ->
            [ meta + metaD + ['subsample_fraction': depth as Float / mean, 'mean_depth': mean, 'depth': depth ], bam, index ]
        }

    ch_input_subsample
        .filter{ meta, _bam, _index -> meta.subsample_fraction >= 1 }
        .subscribe{ meta, _bam, _index ->
            error "Sample ${meta.id} has mean depth of ${meta.mean_depth} and requested depth of ${meta.depth}, resulting in a subsampling fraction of ${meta.subsample_fraction}. As this is >= 1, no subsampling can be performed for this sample."
        }

    // Downsample
    SAMTOOLS_VIEW(
        ch_input_subsample,
        ch_fasta,
        [[], []],
        ch_regions,
        []
    )

    // Aggregate bam and index
    ch_bam_subsampled = SAMTOOLS_VIEW.out.bam
        .mix(SAMTOOLS_VIEW.out.cram, SAMTOOLS_VIEW.out.sam)
        .join(SAMTOOLS_VIEW.out.bai
            .mix(SAMTOOLS_VIEW.out.crai, SAMTOOLS_VIEW.out.csi)
        )
        .map{ meta, bam, index ->
            def keys_to_keep = meta.keySet() - ['subsample_fraction', 'mean_depth']
            [ meta.subMap(keys_to_keep), bam, index ]
        }

    emit:
    bam_subsampled    = ch_bam_subsampled             // channel: [ val(meta), path(bam), path(csi) ]
}
