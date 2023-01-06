//
// Run SAMtools idxstats, then split reads in input bam into one bam per chromosome in idxstats output
//

include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_SPLIT_BY_REGION {

    take:
    ch_bam          // channel: [ val(meta), bam , bai ]
    ch_regions_file // channel: [ regions_file ]

    main:

    ch_versions = Channel.empty()

    // Create channel containing the region names from the bed file.
    ch_regions = ch_regions_file
                    .splitCsv ( header: ['seq_name', 'start', 'stop'], sep:'\t' )
                .map{ stats ->
                    // If the regions file contains just a sequence name provide that
                    if (! stats['start'] ) [ stats['seq_name'] ]
                    // If a specific position is given, use that
                    else if ( ! stats['stop']) [ [ stats['seq_name'], stats['start'] ].join(":") ]
                    // If a region between specific bps is requested use that
                    else [ [ [ stats['seq_name'], stats['start'] ].join(":"), stats['stop'] ].join('-') ]
                }

    // Combine input bam with regions.
    ch_bam
        .combine(ch_regions)
        // Place region into map
        .map{ meta, bam, bai, chrom -> [ meta + [ genomic_region:chrom ], bam, bai ] }
        .set{ ch_bam_for_splitting }

    // The specified region is put into ext.args2 from the meta. See nextflow.config of the subworkflow.
    SAMTOOLS_VIEW(ch_bam_for_splitting, [], [])
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Index the output bams
    SAMTOOLS_INDEX(SAMTOOLS_VIEW.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Emit channel in the same format as was taken in by joining each bam with its bai.
    ch_output = SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_INDEX.out.bai)

    emit:
    bam_bai     = ch_output                         // channel: [ val(meta), bam, bai ]
    versions    = ch_versions                       // channel: [ versions.yml ]
}

