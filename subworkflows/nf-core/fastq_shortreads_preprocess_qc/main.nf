// statistics
include { FASTQ_GENERATE_STATISTICS as PRE_STATS  } from '../fastq_generate_statistics/main'
include { FASTQ_GENERATE_STATISTICS as POST_STATS } from '../fastq_generate_statistics/main'
// preprocessing
include { FASTQ_PREPROCESS_SEQKIT                 } from '../fastq_preprocess_seqkit/main'
// barcoding
include { UMITOOLS_EXTRACT                        } from '../../../modules/nf-core/umitools/extract/main'
// adapter removal and merging
// include { FASTQ_REMOVEADAPTERS_MERGE              } from '../fastq_removeadapters_merge/main'
// complexity filtering
include { PRINSEQPLUSPLUS                         } from '../../../modules/nf-core/prinseqplusplus/main'
// deduplication
include { BBMAP_CLUMPIFY                          } from '../../../modules/nf-core/bbmap/clumpify/main'
// host decontamination
include { FASTQ_DECONTAMINATE_DEACON_HOSTILE      } from '../fastq_decontaminate_deacon_hostile/main'
// final concatenation
include { CAT_FASTQ                               } from '../../../modules/nf-core/cat/fastq/main'

workflow FASTQ_SHORTREADS_PREPROCESS_QC {

    take:
    ch_reads               // channel: [ val(meta), [ fastq ] ]
    // statistics
    skip_fastqc            // boolean
    skip_seqfu_check       // boolean
    skip_seqfu_stats       // boolean
    skip_seqkit_stats      // boolean
    skip_seqtk_comp        // boolean
    // preprocessing
    skip_seqkit_sana_pair  // boolean
    skip_seqkit_seq        // boolean
    skip_seqkit_replace    // boolean
    skip_seqkit_rmdup      // boolean
    // barcoding
    skip_umitools_extract  // boolean
    umi_discard_read       // integer: 0, 1 or 2
    // adapter removal and merging
    // skip_adapterremoval    // boolean
    // complexity filtering
    skip_prinseqplusplus   // boolean
    // deduplication
    skip_bbmap_clumpify    // boolean
    // host decontamination
    skip_decontamination    // boolean
    ch_fasta                // channel: [ val(meta), [ fasta ] ] (optional)
    ch_reference            // channel: [ val(reference_name), path(reference_dir) ] (optional)
    index_name              // val (optional)
    decontaminator          // string (enum): 'hostile' or 'deacon'
    // final concatenation
    skip_cat_fastq          // boolean

    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    // pre-statistics
    PRE_STATS (
        ch_reads,
        skip_fastqc,
        skip_seqfu_check,
        skip_seqfu_stats,
        skip_seqkit_stats,
        skip_seqtk_comp
    )
    ch_versions = ch_versions.mix(PRE_STATS.out.versions)

    // preprocessing
    FASTQ_PREPROCESS_SEQKIT (
        ch_reads,
        skip_seqkit_sana_pair,
        skip_seqkit_seq,
        skip_seqkit_replace,
        skip_seqkit_rmdup
    )
    ch_versions = ch_versions.mix(FASTQ_PREPROCESS_SEQKIT.out.versions)

    // barcoding
    umi_reads = FASTQ_PREPROCESS_SEQKIT.out.reads
    umi_log = channel.empty()
    if (!skip_umitools_extract) {
        UMITOOLS_EXTRACT( ch_reads )
        umi_reads = UMITOOLS_EXTRACT.out.reads
        umi_log = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

        // Discard R1 / R2 if required
        if (umi_discard_read in [1, 2]) {
            UMITOOLS_EXTRACT.out.reads
                .map { meta, reads ->
                    meta.single_end ? [meta, reads] : [meta + ['single_end': true], reads[umi_discard_read % 2]]
                }
                .set { umi_reads }
        }

        ch_reads = umi_reads
    }

    // adapter removal and merging
    // TODO
    // if (!skip_adapterremoval) {

    // }

    // complexity filtering
    // TODO
    // if (!skip_complexity_filtering) {
    //     PRINSEQPLUSPLUS( ... )
    //     ch_versions = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())
    // }

    // deduplication
    // TODO
    // if (!skip_deduplication) {
    //     BBMAP_CLUMPIFY( ... )
    //     ch_versions = ch_versions.mix(BBMAP_CLUMPIFY.out.versions.first())
    // }

    // host decontamination
    if (!skip_decontamination) {
        FASTQ_DECONTAMINATE_DEACON_HOSTILE (
            ch_reads,
            ch_fasta,
            ch_reference,
            index_name,
            decontaminator
        )
        ch_versions = ch_versions.mix(FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.versions)

        ch_reads = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.fastq_filtered
    }


    // final concatenation
    // TODO
    // if (!skip_final_concatenation) {
    //     CAT_FASTQ( ... )
    //     ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
    // }

    // post-statistics
    POST_STATS (
        ch_reads,
        skip_fastqc,
        skip_seqfu_check,
        skip_seqfu_stats,
        skip_seqkit_stats,
        skip_seqtk_comp
    )
    ch_versions = ch_versions.mix(POST_STATS.out.versions)

    emit:
    reads         = ch_reads          // channel: [ val(meta), [ fastq ] ]
    versions      = ch_versions       // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
