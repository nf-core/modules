//
// Estimate telomere length and content from aligned reads
//
// Length estimation: telseq (short reads) or telogator2 (long reads)
// Content profiling: telomerehunter (optional, independent of length estimation)
// Auto-indexes BAMs/CRAMs and FASTAs when indices not provided
//

include { SAMTOOLS_INDEX                     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX                     } from '../../../modules/nf-core/samtools/faidx/main'
include { TELSEQ                             } from '../../../modules/nf-core/telseq/main'
include { TELOMEREHUNTER                     } from '../../../modules/nf-core/telomerehunter/main'
include { TELOGATOR2                         } from '../../../modules/nf-core/telogator2/main'
include { CUSTOM_SUMMARISETELOMEREESTIMATION } from '../../../modules/nf-core/custom/summarisetelomereestimation/main'

def classifyReadType(data_type) {
    if (data_type in ['wgs', 'wes', 'short_read', 'telseq'])       return 'short_read'
    if (data_type in ['long_read', 'ont', 'pacbio', 'telogator2']) return 'long_read'
    return null
}

workflow BAM_TELOMERE_ESTIMATION {

    take:
    ch_reads            // channel: [ val(meta), val(data_type), path(reads), path(reads_index) ]
    ch_control          // channel: [ val(meta), path(control_bam), path(control_bai) ]
    ch_fasta            // channel: [ val(meta2), path(fasta), path(fai) ]
    ch_bed              // channel: [ val(meta4), path(bed) ]
    run_telomerehunter  // val(boolean): enable telomerehunter content profiling
    length_estimator    // val(string):  'telseq'|'telogator2'|null - global override

    main:

    //
    // AUTO-INDEX: FASTA
    //
    ch_fasta_branched = ch_fasta.branch { _meta, fasta, fai ->
        indexed:     fai
        needs_index: fasta
        empty:       true
    }

    SAMTOOLS_FAIDX(
        ch_fasta_branched.needs_index,
        false
    )

    ch_fasta_ready = ch_fasta_branched.indexed
        .mix(ch_fasta_branched.empty)
        .mix(
            ch_fasta_branched.needs_index
                .join(SAMTOOLS_FAIDX.out.fai)
                .map { meta, fasta, _old_fai, fai -> [ meta, fasta, fai ] }
        )

    //
    // AUTO-INDEX: BAM/CRAM
    //
    ch_reads_idx_branched = ch_reads.branch { _meta, _dt, _reads, index ->
        indexed:   index
        unindexed: true
    }

    SAMTOOLS_INDEX(
        ch_reads_idx_branched.unindexed.map { meta, _dt, reads, _index -> [ meta, reads ] }
    )

    ch_reads_with_index = ch_reads_idx_branched.indexed
        .mix(
            ch_reads_idx_branched.unindexed
                .join(SAMTOOLS_INDEX.out.index)
                .map { meta, dt, reads, _old_idx, idx -> [ meta, dt, reads, idx ] }
        )

    //
    // Route by data_type
    //
    // Global length_estimator overrides per-sample data_type when set.
    // Samples without a resolvable data_type will error.
    //
    ch_reads_categorised = ch_reads_with_index.map { meta, data_type, reads, index ->
        def effective_type = length_estimator ?: data_type
        def category = classifyReadType(effective_type)
        if (category == null) {
            error("Unrecognised data_type '${effective_type}' for sample ${meta.id}. " +
                  "Valid values: wgs, wes, short_read, telseq, long_read, ont, pacbio, telogator2.")
        }
        [ meta + [read_category: category], reads, index ]
    }

    //
    // TELSEQ: short-read telomere length
    //
    TELSEQ(
        ch_reads_categorised.filter { meta, _reads, _index -> meta.read_category == 'short_read' },
        ch_fasta_ready,
        ch_bed
    )

    //
    // TELOGATOR2: long-read telomere length
    //
    TELOGATOR2(
        ch_reads_categorised.filter { meta, _reads, _index -> meta.read_category == 'long_read' },
        ch_fasta_ready
    )

    //
    // Collect length estimation TSVs
    //
    ch_length_tsv = TELSEQ.out.output
        .mix(TELOGATOR2.out.tlens)

    //
    // TELOMEREHUNTER: optional content profiling
    //
    // Pair each sample with its matched control (if any), then gate on run_telomerehunter
    ch_th_input = ch_reads_categorised
        .join(ch_control, remainder: true)
        .map { items ->
            def meta = items[0]
            def bam  = items[1]
            def bai  = items[2]
            def cbam = items.size() > 3 ? items[3] : []
            def cbai = items.size() > 4 ? items[4] : []
            [ meta, bam, bai, cbam ?: [], cbai ?: [] ]
        }
        .filter { _meta, _bam, _bai, _cbam, _cbai -> run_telomerehunter as boolean }

    TELOMEREHUNTER(ch_th_input)

    ch_content_tsv          = TELOMEREHUNTER.out.summary
    ch_telomerehunter_tumor = TELOMEREHUNTER.out.tumor
    ch_telomerehunter_ctrl  = TELOMEREHUNTER.out.control

    //
    // Summarise: left-join length results with optional content results per sample
    //
    ch_summary_input = ch_length_tsv
        .join(ch_content_tsv, remainder: true)
        .map { items ->
            def meta        = items[0]
            def length_tsv  = items[1]
            def content_tsv = items.size() > 2 ? items[2] : []
            def tool = meta.read_category == 'short_read' ? 'telseq' : 'telogator2'
            [ meta, length_tsv, content_tsv ?: [], tool ]
        }

    CUSTOM_SUMMARISETELOMEREESTIMATION(
        ch_summary_input
    )

    emit:
    length_tsv             = ch_length_tsv                                      // channel: [ val(meta), path(tsv) ]
    content_tsv            = ch_content_tsv                                     // channel: [ val(meta), path(tsv) ]
    summary                = CUSTOM_SUMMARISETELOMEREESTIMATION.out.summary     // channel: [ val(meta), path(tsv) ]
    telomerehunter_tumor   = ch_telomerehunter_tumor                            // channel: [ val(meta), path(dir) ]
    telomerehunter_control = ch_telomerehunter_ctrl                             // channel: [ val(meta), path(dir) ]
    telogator2_plots       = TELOGATOR2.out.plots                               // channel: [ val(meta), path(*.png) ]
    telogator2_stats       = TELOGATOR2.out.stats                               // channel: [ val(meta), path(stats.txt) ]
}
