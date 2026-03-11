//
// Run SAMtools idxstats, then split reads in input bam into one bam per chromosome in idxstats output
//

include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_SPLIT_BY_REGION {

    take:
    ch_bam // channel: [ val(meta), path(bam), path(bai), path(regions_file) ]

    main:

    ch_versions = channel.empty()

    //
    // Create channel containing the region names from the bed file.
    //

    ch_regions = ch_bam
        .map{
            meta, _bam, _bai, regions_file ->
            [ meta, regions_file ]
        }
        .splitCsv ( header: ['seq_name', 'start', 'stop'], sep:'\t', elem: 1)
        .map{ meta, stats ->
            // If the regions file contains just a sequence name provide that
            if (!stats['start']) {
                return [meta, stats['seq_name']]
            }

            // If a specific position is given, use that
            if (!stats['stop']) {
                return [meta, [stats['seq_name'], stats['start']].join(":")]
            }

            // If a region between specific bps is requested use that
            def chrom = [[stats['seq_name'], stats['start']].join(":"), stats['stop']].join('-')
            return [meta, chrom]
        }

    //
    // Combine input bam with regions.
    //

    ch_bam_for_splitting = ch_bam
        .combine(ch_regions, by: 0)
        // Place region into meta map
        .map{ meta, bam, bai, _region_file, chrom -> [ meta + [ genomic_region:chrom ], bam, bai ] }

    // The specified region is put into ext.args2 from the meta. See nextflow.config of the subworkflow.
    SAMTOOLS_VIEW(ch_bam_for_splitting, [[],[],[]], [], [])

    //
    // Index the output bams
    //

    SAMTOOLS_INDEX(SAMTOOLS_VIEW.out.bam)

    //
    // Emit channel in the same format as was taken in by joining each bam with its bai.
    //

    ch_output = SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_INDEX.out.bai)

    emit:
    bam_bai     = ch_output                         // channel: [ val(meta), path(bam), path(bai) ]
    versions    = ch_versions                       // channel: [ path(versions.yml) ]
}
