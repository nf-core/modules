// statistics
include { FASTQ_QC_STATS as PRE_STATS        } from '../fastq_qc_stats/main'
include { FASTQ_QC_STATS as POST_STATS       } from '../fastq_qc_stats/main'
// preprocessing
include { FASTQ_PREPROCESS_SEQKIT            } from '../fastq_preprocess_seqkit/main'
// barcoding
include { UMITOOLS_EXTRACT                   } from '../../../modules/nf-core/umitools/extract/main'
// adapter removal and merging
include { FASTQ_REMOVEADAPTERS_MERGE         } from '../fastq_removeadapters_merge/main'
// complexity filtering
include { PRINSEQPLUSPLUS                    } from '../../../modules/nf-core/prinseqplusplus/main'
// deduplication
include { BBMAP_CLUMPIFY                     } from '../../../modules/nf-core/bbmap/clumpify/main'
// host decontamination
include { FASTQ_DECONTAMINATE_DEACON_HOSTILE } from '../fastq_decontaminate_deacon_hostile/main'
// final concatenation
include { CAT_FASTQ                          } from '../../../modules/nf-core/cat/fastq/main'

workflow FASTQ_SHORTREADS_PREPROCESS_QC {

    take:
    ch_reads                       // channel: [ val(meta), [ fastq ] ]
    // statistics
    skip_fastqc                    // boolean
    skip_seqfu_check               // boolean
    skip_seqfu_stats               // boolean
    skip_seqkit_stats              // boolean
    skip_seqtk_comp                // boolean
    // preprocessing
    skip_seqkit_sana_pair          // boolean
    skip_seqkit_seq                // boolean
    skip_seqkit_replace            // boolean
    skip_seqkit_rmdup              // boolean
    // barcoding
    skip_umitools_extract          // boolean
    val_umi_discard_read           // integer: 0, 1 or 2
    // adapter removal and merging
    skip_adapterremoval            // boolean
    val_adapter_tool               // string:  [mandatory] tool_name // choose from: ["trimmomatic", "cutadapt", "trimgalore", "bbduk", "leehom", "fastp", "adapterremoval"]
    ch_custom_adapters_file        // channel: [optional]  {fasta,txt} // fasta, for bbduk or fastp, or txt, for adapterremoval
    val_save_merged                // boolean: [mandatory] if true, will return the merged reads instead, for fastp and adapterremoval
    val_fastp_discard_trimmed_pass // boolean: [mandatory] // only for fastp
    val_fastp_save_trimmed_fail    // boolean: [mandatory] // only for fastp
    // complexity filtering
    skip_complexity_filtering      // boolean
    // deduplication
    skip_deduplication             // boolean
    // host decontamination
    skip_decontamination           // boolean
    ch_fasta                       // channel: [ val(meta), [ fasta ] ] (optional)
    ch_reference                   // channel: [ val(reference_name), path(reference_dir) ] (optional)
    val_index_name                 // val (optional)
    val_decontaminator             // string (enum): 'hostile' or 'deacon'
    // final concatenation
    skip_final_concatenation       // boolean

    main:

    ch_versions                       = channel.empty()
    ch_multiqc_files                  = channel.empty()
    ch_pre_stats_fastqc_html          = channel.empty()
    ch_pre_stats_fastqc_zip           = channel.empty()
    ch_pre_stats_seqfu_check          = channel.empty()
    ch_pre_stats_seqfu_stats          = channel.empty()
    ch_pre_stats_seqkit_stats         = channel.empty()
    ch_pre_stats_seqtk_stats          = channel.empty()
    ch_post_stats_fastqc_html         = channel.empty()
    ch_post_stats_fastqc_zip          = channel.empty()
    ch_post_stats_seqfu_check         = channel.empty()
    ch_post_stats_seqfu_stats         = channel.empty()
    ch_post_stats_seqkit_stats        = channel.empty()
    ch_post_stats_seqtk_stats         = channel.empty()
    ch_umi_log                        = channel.empty()
    ch_adapterremoval_discarded_reads = channel.empty()
    ch_adapterremoval_logfile         = channel.empty()
    ch_adapterremoval_report          = channel.empty()
    ch_prinseq_log                    = channel.empty()
    ch_clumpify_log                   = channel.empty()
    ch_hostile_reference              = channel.empty()
    ch_hostile_json                   = channel.empty()
    ch_deacon_index                   = channel.empty()
    ch_deacon_summary                 = channel.empty()

    // pre-statistics
    PRE_STATS (
        ch_reads,
        skip_fastqc,
        skip_seqfu_check,
        skip_seqfu_stats,
        skip_seqkit_stats,
        skip_seqtk_comp
    )
    ch_pre_stats_fastqc_html  = PRE_STATS.out.fastqc_html
    ch_pre_stats_fastqc_zip   = PRE_STATS.out.fastqc_zip
    ch_pre_stats_seqfu_check  = PRE_STATS.out.seqfu_check
    ch_pre_stats_seqfu_stats  = PRE_STATS.out.seqfu_stats
    ch_pre_stats_seqkit_stats = PRE_STATS.out.seqkit_stats
    ch_pre_stats_seqtk_stats  = PRE_STATS.out.seqtk_stats
    ch_multiqc_files          = ch_multiqc_files.mix(PRE_STATS.out.seqfu_multiqc)
    ch_versions               = ch_versions.mix(PRE_STATS.out.versions)

    // preprocessing
    FASTQ_PREPROCESS_SEQKIT (
        ch_reads,
        skip_seqkit_sana_pair,
        skip_seqkit_seq,
        skip_seqkit_replace,
        skip_seqkit_rmdup
    )
    ch_reads    = FASTQ_PREPROCESS_SEQKIT.out.reads
    ch_versions = ch_versions.mix(FASTQ_PREPROCESS_SEQKIT.out.versions)

    // barcoding
    if (!skip_umitools_extract) {
        UMITOOLS_EXTRACT( ch_reads )
        ch_umi_reads = UMITOOLS_EXTRACT.out.reads
        ch_umi_log = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

        // Discard R1 / R2 if required
        if (val_umi_discard_read in [1, 2]) {
            ch_umi_reads = UMITOOLS_EXTRACT.out.reads
                .map { meta, reads ->
                    meta.single_end ? [meta, reads] : [meta + ['single_end': true], reads[val_umi_discard_read % 2]]
                }
        }

        ch_reads = ch_umi_reads
    }

    // adapter removal and merging
    if (!skip_adapterremoval) {
        FASTQ_REMOVEADAPTERS_MERGE (
            ch_reads,
            val_adapter_tool,
            ch_custom_adapters_file,
            val_save_merged,
            val_fastp_discard_trimmed_pass,
            val_fastp_save_trimmed_fail
        )
        ch_adapterremoval_discarded_reads = FASTQ_REMOVEADAPTERS_MERGE.out.discarded_reads
        ch_adapterremoval_logfile         = FASTQ_REMOVEADAPTERS_MERGE.out.ch_log
        ch_adapterremoval_report          = FASTQ_REMOVEADAPTERS_MERGE.out.ch_report
        ch_reads                          = FASTQ_REMOVEADAPTERS_MERGE.out.processed_reads
        ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_REMOVEADAPTERS_MERGE.out.multiqc_files)
        ch_versions                       = ch_versions.mix(FASTQ_REMOVEADAPTERS_MERGE.out.versions)
    }

    // complexity filtering
    if (!skip_complexity_filtering) {
        PRINSEQPLUSPLUS( ch_reads )
        ch_reads       = PRINSEQPLUSPLUS.out.good_reads
        ch_prinseq_log = PRINSEQPLUSPLUS.out.log
        ch_versions    = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())
    }

    // deduplication
    if (!skip_deduplication) {
        BBMAP_CLUMPIFY( ch_reads )
        ch_reads        = BBMAP_CLUMPIFY.out.reads
        ch_clumpify_log = BBMAP_CLUMPIFY.out.log
        ch_versions     = ch_versions.mix(BBMAP_CLUMPIFY.out.versions.first())
    }

    // host decontamination
    if (!skip_decontamination) {
        FASTQ_DECONTAMINATE_DEACON_HOSTILE (
            ch_reads,
            ch_fasta,
            ch_reference,
            val_index_name,
            val_decontaminator
        )
        ch_reads             = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.fastq_filtered
        ch_hostile_reference = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.reference
        ch_hostile_json      = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.json
        ch_deacon_index      = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.index
        ch_deacon_summary    = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.summary
        ch_versions          = ch_versions.mix(FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.versions)
    }


    // final concatenation
    if (!skip_final_concatenation) {
        // CAT_FASTQ ( ch_reads.map { meta, reads -> [meta, reads.flatten()] } ) // TODO test more cases
        CAT_FASTQ ( ch_reads )
        ch_reads = CAT_FASTQ.out.reads
    }

    // post-statistics
    POST_STATS (
        ch_reads,
        skip_fastqc,
        skip_seqfu_check,
        skip_seqfu_stats,
        skip_seqkit_stats,
        skip_seqtk_comp
    )
    ch_post_stats_fastqc_html  = POST_STATS.out.fastqc_html
    ch_post_stats_fastqc_zip   = POST_STATS.out.fastqc_zip
    ch_post_stats_seqfu_check  = POST_STATS.out.seqfu_check
    ch_post_stats_seqfu_stats  = POST_STATS.out.seqfu_stats
    ch_post_stats_seqkit_stats = POST_STATS.out.seqkit_stats
    ch_post_stats_seqtk_stats  = POST_STATS.out.seqtk_stats
    ch_multiqc_files           = ch_multiqc_files.mix(POST_STATS.out.seqfu_multiqc)
    ch_versions                = ch_versions.mix(POST_STATS.out.versions)

    emit:
    reads = ch_reads // channel: [ val(meta), [ fastq ] ]

    // statistics
    pre_stats_fastqc_html    = ch_pre_stats_fastqc_html
    pre_stats_fastqc_zip     = ch_pre_stats_fastqc_zip
    pre_stats_seqfu_check    = ch_pre_stats_seqfu_check
    pre_stats_seqfu_stats    = ch_pre_stats_seqfu_stats
    pre_stats_seqkit_stats   = ch_pre_stats_seqkit_stats
    pre_stats_seqtk_stats    = ch_pre_stats_seqtk_stats
    post_stats_fastqc_html   = ch_post_stats_fastqc_html
    post_stats_fastqc_zip    = ch_post_stats_fastqc_zip
    post_stats_seqfu_check   = ch_post_stats_seqfu_check
    post_stats_seqfu_stats   = ch_post_stats_seqfu_stats
    post_stats_seqkit_stats  = ch_post_stats_seqkit_stats
    post_stats_seqtk_stats   = ch_post_stats_seqtk_stats

    // barcoding
    umi_log = ch_umi_log

    // adapter removal and merging
    adapterremoval_discarded_reads = ch_adapterremoval_discarded_reads
    adapterremoval_logfile         = ch_adapterremoval_logfile
    adapterremoval_report          = ch_adapterremoval_report

    // complexity filtering
    prinseq_log = ch_prinseq_log

    // deduplication
    clumpify_log = ch_clumpify_log

    // host decontamination
    hostile_reference = ch_hostile_reference
    hostile_json      = ch_hostile_json
    deacon_index      = ch_deacon_index
    deacon_summary    = ch_deacon_summary

    versions      = ch_versions      // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
