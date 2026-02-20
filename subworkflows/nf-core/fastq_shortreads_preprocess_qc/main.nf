// statistics
include { FASTQ_QC_STATS as PRE_STATS        } from '../fastq_qc_stats'
include { FASTQ_QC_STATS as POST_STATS       } from '../fastq_qc_stats'
// preprocessing
include { FASTQ_PREPROCESS_SEQKIT            } from '../fastq_preprocess_seqkit'
// barcoding
include { UMITOOLS_EXTRACT                   } from '../../../modules/nf-core/umitools/extract'
// adapter removal and merging
include { FASTQ_REMOVEADAPTERS_MERGE         } from '../fastq_removeadapters_merge'
// complexity filtering
include { FASTQ_COMPLEXITY_FILTER            } from '../fastq_complexity_filter'
// deduplication
include { BBMAP_CLUMPIFY                     } from '../../../modules/nf-core/bbmap/clumpify'
// host decontamination
include { FASTQ_DECONTAMINATE_DEACON_HOSTILE } from '../fastq_decontaminate_deacon_hostile'
// final concatenation
include { CAT_FASTQ                          } from '../../../modules/nf-core/cat/fastq'

workflow FASTQ_SHORTREADS_PREPROCESS_QC {
    take:
    ch_reads // channel: [ val(meta), [ fastq ] ]
    skip_fastqc // boolean
    skip_seqfu_check // boolean
    skip_seqfu_stats // boolean
    skip_seqkit_stats // boolean
    skip_seqtk_comp // boolean
    skip_seqkit_sana_pair // boolean
    skip_seqkit_seq // boolean
    skip_seqkit_replace // boolean
    skip_seqkit_rmdup // boolean
    skip_umitools_extract // boolean
    val_umi_discard_read // integer: 0, 1 or 2
    skip_adapterremoval // boolean
    val_adapter_tool // string:  [mandatory] tool_name // choose from: ["trimmomatic", "cutadapt", "trimgalore", "bbduk", "leehom", "fastp", "adapterremoval"]
    ch_custom_adapters_file // channel: [optional]  [ {fasta,txt} ] // fasta, for bbduk or fastp, or txt, for adapterremoval
    val_save_merged // boolean: [mandatory] if true, will return the merged reads instead, for fastp and adapterremoval
    val_fastp_discard_trimmed_pass // boolean: [mandatory] // only for fastp
    val_fastp_save_trimmed_fail // boolean: [mandatory] // only for fastp
    skip_complexity_filtering // boolean
    val_complexity_filter_tool // string:  [mandatory] tool_name // choose from: ["prinseqplusplus", "bbduk", "fastp"]
    skip_deduplication // boolean
    skip_decontamination // boolean
    ch_decontamination_fasta // channel: [ val(meta), [ fasta ] ] (optional)
    ch_decontamination_reference // channel: [ val(reference_name), path(reference_dir) ] (optional)
    val_decontamination_index_name // val (optional)
    val_decontamination_tool // string (enum): 'hostile' or 'deacon'
    skip_final_concatenation // boolean

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_umi_log = channel.empty()
    ch_adapterremoval_discarded_reads = channel.empty()
    ch_adapterremoval_logfile = channel.empty()
    ch_adapterremoval_report = channel.empty()
    ch_complexity_filter_log = channel.empty()
    ch_complexity_filter_report = channel.empty()
    ch_clumpify_log = channel.empty()
    ch_hostile_reference = channel.empty()
    ch_hostile_json = channel.empty()
    ch_deacon_index = channel.empty()
    ch_deacon_summary = channel.empty()

    // pre-statistics
    PRE_STATS(
        ch_reads,
        skip_fastqc,
        skip_seqfu_check,
        skip_seqfu_stats,
        skip_seqkit_stats,
        skip_seqtk_comp,
    )
    ch_pre_stats_fastqc_html = PRE_STATS.out.fastqc_html
    ch_pre_stats_fastqc_zip = PRE_STATS.out.fastqc_zip
    ch_pre_stats_seqfu_check = PRE_STATS.out.seqfu_check
    ch_pre_stats_seqfu_stats = PRE_STATS.out.seqfu_stats
    ch_pre_stats_seqkit_stats = PRE_STATS.out.seqkit_stats
    ch_pre_stats_seqtk_stats = PRE_STATS.out.seqtk_stats
    ch_multiqc_files = ch_multiqc_files.mix(PRE_STATS.out.seqfu_multiqc)

    // preprocessing
    FASTQ_PREPROCESS_SEQKIT(
        ch_reads,
        skip_seqkit_sana_pair,
        skip_seqkit_seq,
        skip_seqkit_replace,
        skip_seqkit_rmdup,
    )
    ch_reads = FASTQ_PREPROCESS_SEQKIT.out.reads
    ch_versions = ch_versions.mix(FASTQ_PREPROCESS_SEQKIT.out.versions)

    // barcoding
    if (!skip_umitools_extract) {
        UMITOOLS_EXTRACT(ch_reads)
        ch_umi_reads = UMITOOLS_EXTRACT.out.reads
        ch_umi_log = UMITOOLS_EXTRACT.out.log

        // Discard R1 / R2 if required
        if (val_umi_discard_read in [1, 2]) {
            ch_umi_reads = UMITOOLS_EXTRACT.out.reads.map { meta, reads ->
                meta.single_end ? [meta, reads] : [meta + ['single_end': true], reads[val_umi_discard_read % 2]]
            }
        }

        ch_reads = ch_umi_reads
    }

    // adapter removal and merging
    if (!skip_adapterremoval) {
        FASTQ_REMOVEADAPTERS_MERGE(
            ch_reads,
            val_adapter_tool,
            ch_custom_adapters_file,
            val_save_merged,
            val_fastp_discard_trimmed_pass,
            val_fastp_save_trimmed_fail,
        )
        ch_adapterremoval_discarded_reads = FASTQ_REMOVEADAPTERS_MERGE.out.discarded_reads
        ch_adapterremoval_logfile = FASTQ_REMOVEADAPTERS_MERGE.out.logfile
        ch_adapterremoval_report = FASTQ_REMOVEADAPTERS_MERGE.out.report
        ch_reads = FASTQ_REMOVEADAPTERS_MERGE.out.processed_reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_REMOVEADAPTERS_MERGE.out.multiqc_files)
        ch_versions = ch_versions.mix(FASTQ_REMOVEADAPTERS_MERGE.out.versions)
    }

    // complexity filtering
    if (!skip_complexity_filtering) {
        FASTQ_COMPLEXITY_FILTER(ch_reads, val_complexity_filter_tool)
        ch_reads = FASTQ_COMPLEXITY_FILTER.out.filtered_reads
        ch_complexity_filter_log = FASTQ_COMPLEXITY_FILTER.out.logfile
        ch_complexity_filter_report = FASTQ_COMPLEXITY_FILTER.out.report
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_COMPLEXITY_FILTER.out.multiqc_files)
        ch_versions = ch_versions.mix(FASTQ_COMPLEXITY_FILTER.out.versions)
    }

    // deduplication
    if (!skip_deduplication) {
        BBMAP_CLUMPIFY(ch_reads)
        ch_reads = BBMAP_CLUMPIFY.out.reads
        ch_clumpify_log = BBMAP_CLUMPIFY.out.log
        ch_versions = ch_versions.mix(BBMAP_CLUMPIFY.out.versions.first())
    }

    // host decontamination
    if (!skip_decontamination) {
        FASTQ_DECONTAMINATE_DEACON_HOSTILE(
            ch_reads,
            ch_decontamination_fasta,
            ch_decontamination_reference,
            val_decontamination_index_name,
            val_decontamination_tool,
        )
        ch_reads = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.fastq_filtered
        ch_hostile_reference = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.reference
        ch_hostile_json = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.json
        ch_deacon_index = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.index
        ch_deacon_summary = FASTQ_DECONTAMINATE_DEACON_HOSTILE.out.summary
    }


    // final concatenation
    if (!skip_final_concatenation) {
        ch_reads_for_cat_branch = ch_reads
            .groupTuple()
            .map { meta, reads ->
                [meta, reads.flatten()]
            }
            .branch { meta, reads ->
                cat: (meta.single_end && reads.size() > 1) || (!meta.single_end && reads.size() > 2)
                skip: true
            }

        CAT_FASTQ(ch_reads_for_cat_branch.cat)

        ch_reads = CAT_FASTQ.out.reads
            .mix(ch_reads_for_cat_branch.skip)
            .map { meta, reads ->
                def new_reads = meta.single_end ? reads[0] : reads.flatten()
                [meta, new_reads]
            }
    }

    // post-statistics
    POST_STATS(
        ch_reads,
        skip_fastqc,
        skip_seqfu_check,
        skip_seqfu_stats,
        skip_seqkit_stats,
        skip_seqtk_comp,
    )
    ch_post_stats_fastqc_html = POST_STATS.out.fastqc_html
    ch_post_stats_fastqc_zip = POST_STATS.out.fastqc_zip
    ch_post_stats_seqfu_check = POST_STATS.out.seqfu_check
    ch_post_stats_seqfu_stats = POST_STATS.out.seqfu_stats
    ch_post_stats_seqkit_stats = POST_STATS.out.seqkit_stats
    ch_post_stats_seqtk_stats = POST_STATS.out.seqtk_stats
    ch_multiqc_files = ch_multiqc_files.mix(POST_STATS.out.seqfu_multiqc)

    emit:
    reads                          = ch_reads // channel: [ val(meta), [ fastq ] ]
    pre_stats_fastqc_html          = ch_pre_stats_fastqc_html
    pre_stats_fastqc_zip           = ch_pre_stats_fastqc_zip
    pre_stats_seqfu_check          = ch_pre_stats_seqfu_check
    pre_stats_seqfu_stats          = ch_pre_stats_seqfu_stats
    pre_stats_seqkit_stats         = ch_pre_stats_seqkit_stats
    pre_stats_seqtk_stats          = ch_pre_stats_seqtk_stats
    post_stats_fastqc_html         = ch_post_stats_fastqc_html
    post_stats_fastqc_zip          = ch_post_stats_fastqc_zip
    post_stats_seqfu_check         = ch_post_stats_seqfu_check
    post_stats_seqfu_stats         = ch_post_stats_seqfu_stats
    post_stats_seqkit_stats        = ch_post_stats_seqkit_stats
    post_stats_seqtk_stats         = ch_post_stats_seqtk_stats
    umi_log                        = ch_umi_log
    adapterremoval_discarded_reads = ch_adapterremoval_discarded_reads
    adapterremoval_logfile         = ch_adapterremoval_logfile
    adapterremoval_report          = ch_adapterremoval_report
    complexity_filter_log          = ch_complexity_filter_log
    complexity_filter_report       = ch_complexity_filter_report
    clumpify_log                   = ch_clumpify_log
    hostile_reference              = ch_hostile_reference
    hostile_json                   = ch_hostile_json
    deacon_index                   = ch_deacon_index
    deacon_summary                 = ch_deacon_summary
    multiqc_files                  = ch_multiqc_files
    versions                       = ch_versions // channel: [ versions.yml ]
}
