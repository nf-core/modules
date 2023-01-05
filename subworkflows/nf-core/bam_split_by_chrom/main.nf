//
// Run SAMtools idxstats, then split reads in input bam into one bam per chromosome in idxstats output
//

include { SAMTOOLS_IDXSTATS      } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_VIEW          } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_SPLIT_BY_CHROM {

    take:
    ch_bam // channel: [ val(meta), bam , bai ]

    main:

    ch_versions = Channel.empty()

    // Run idxstats to get names for chromosomes but also umber of mapped reads.
    SAMTOOLS_IDXSTATS(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    // Create channel containing the sequence names that contain any mapped reads in them for each meta entry. [ [meta], seq_name ]
    ch_chroms_with_reads = SAMTOOLS_IDXSTATS.out.idxstats
                                .splitCsv ( header: ['seq_name', 'seq_length', 'n_mapped_reads', 'n_unmapped_reads'], sep:'\t' )
                                // Only keep chromosomes that have some data to avoid creating empty bams
                                //      that will likely fail in some downstream processes
                                .filter{ it[1]['n_mapped_reads'].toInteger() > 0 }
                                .map{
                                    meta, stats ->
                                    [ meta, stats['seq_name'] ]
                                }

    // Create a channel for each input bam that has an additional field iterating over all chromosomes containing reads
    ch_bam
        .combine(ch_chroms_with_reads, by:0)
        // Place chromosome name into map
        .map{
            meta, bam, bai, chrom ->
            clone = meta.clone()
            clone['seq_name'] = chrom

            [ clone, bam, bai ]
        }
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
    idxstats    = SAMTOOLS_IDXSTATS.out.idxstats    // channel: [ val(meta), idxstats ]
    versions    = ch_versions                       // channel: [ versions.yml ]
}

