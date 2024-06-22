//
// Alignment with mapAD and sort
//

include { MAPAD_MAP               } from '../../../modules/nf-core/mapad/map/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_MAPAD {

    take:
    ch_reads // channel (mandatory): [ val(meta), path(reads) ]. Important: meta REQUIRES single_end` entry!
    ch_index // channel (mandatory): [ val(meta), path(index) ]
    ch_fasta // channel (optional) : [ val(meta3), path(fasta) ]
    val_mismatch_parameter
    val_double_stranded_library
    val_five_prime_overhang
    val_three_prime_overhang
    val_deam_rate_double_stranded
    val_deam_rate_single_stranded
    val_indel_rate


    main:

    ch_versions = Channel.empty()

    // WARNING: You must specify in your prefix `meta.id_index` in your `modules.conf`
    // to ensure that you do not overwrite multiple BAM files from one sample mapped
    // against multiple references. This meta map is added by the subworkflow but can be removed
    // after if necessary.

    // Ensure when multiple references are provided, that reference/read combinations
    // are correctly associated throughout the subworkflow by copying the sample
    // specific metadata to the index on each combination

    ch_prepped_input = ch_reads
                        .combine(ch_index)
                        .map{
                            meta, reads, meta_index, index ->

                                // Create a combined id that includes the ids of the reads and the index used.
                                // Also keep the id of the index with a new name to avoid name collisions.
                                def key_read_ref = meta.id + "_" + meta_index.id
                                def id_index = meta_index.id

                            [ meta + [key_read_ref: key_read_ref] + [id_index: id_index], reads, meta_index + [key_read_ref: key_read_ref]  + [id_index: id_index], index  ]
                        }

    // Drop the index_meta, as the id of the index is now kept within the read meta.
    ch_preppedinput_for_mapad = ch_prepped_input
                        .multiMap {
                            meta, reads, meta_index, index ->
                                reads: [ meta, reads ]
                                index: [ meta, index ]
                        }

    // Alignment
    MAPAD_MAP ( ch_preppedinput_for_mapad.reads, ch_preppedinput_for_mapad.index, val_mismatch_parameter, val_double_stranded_library, val_five_prime_overhang, val_three_prime_overhang, val_deam_rate_double_stranded, val_deam_rate_single_stranded, val_indel_rate )
    ch_versions = ch_versions.mix( MAPAD_MAP.out.versions.first() )

    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    BAM_SORT_STATS_SAMTOOLS ( MAPAD_MAP.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam_unsorted = MAPAD_MAP.out.bam                    // channel: [ val(meta), path(bam) ]

    bam          = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bai          = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    csi          = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), path(csi) ]
    stats        = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat     = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats     = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                              // channel: [ path(versions.yml) ]
}
