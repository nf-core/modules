//
// Estimate telomere length and content from aligned reads
//
// Length estimation: telseq (short reads) or telogator2 (long reads)
// Content profiling: telomerehunter (optional, independent of length estimation)
//

include { TELSEQ                             } from '../../../modules/nf-core/telseq/main'
include { TELOMEREHUNTER                     } from '../../../modules/nf-core/telomerehunter/main'
include { TELOGATOR2                         } from '../../../modules/nf-core/telogator2/main'
include { CUSTOM_SUMMARISETELOMEREESTIMATION } from '../../../modules/nf-core/custom/summarisetelomereestimation/main'

workflow BAM_TELOMERE_ESTIMATION {

    take:
    ch_reads            // channel: [ val(meta), val(data_type), path(reads), path(reads_index), path(bed) ]
                        //   data_type: 'short' or 'long'; bed: optional exome BED for telseq ([] if not needed)
    ch_control          // channel: [ val(meta), path(control_bam), path(control_bai) ]
    ch_fasta            // channel: [ val(meta2), path(fasta), path(fai) ]
    ch_cytoband         // channel: [ path(cytoband) ] - optional cytoband file for telomerehunter
    run_telomerehunter  // val(boolean): enable telomerehunter content profiling
    length_estimator    // val(string):  'short'|'long'|null - global override

    main:

    //
    // Route by data_type ('short' -> telseq, 'long' -> telogator2)
    //
    ch_reads_validated = ch_reads.map { meta, data_type, reads, index, bed ->
        def category = length_estimator ?: data_type
        if (!(category in ['short', 'long'])) {
            error("Invalid data_type '${category}' for sample ${meta.id}. Must be 'short' or 'long'.")
        }
        [ meta, category, reads, index, bed ]
    }

    ch_short_reads = ch_reads_validated.filter { _meta, dt, _reads, _index, _bed -> dt == 'short' }
    ch_long_reads  = ch_reads_validated.filter { _meta, dt, _reads, _index, _bed -> dt == 'long' }

    //
    // TELSEQ: short-read telomere length
    //
    TELSEQ(
        ch_short_reads.map { meta, _dt, reads, index, bed -> [ meta, reads, index, bed ] },
        ch_fasta
    )

    //
    // TELOGATOR2: long-read telomere length
    //
    TELOGATOR2(
        ch_long_reads.map { meta, _dt, reads, index, _bed -> [ meta, reads, index ] },
        ch_fasta
    )

    //
    // Collect length estimation TSVs, tagging with the tool that produced them
    //
    ch_length_tsv = TELSEQ.out.output
        .mix(TELOGATOR2.out.tlens)

    ch_length_with_tool = TELSEQ.out.output.map { meta, tsv -> [ meta, tsv, 'telseq' ] }
        .mix(TELOGATOR2.out.tlens.map { meta, tsv -> [ meta, tsv, 'telogator2' ] })

    //
    // TELOMEREHUNTER: optional content profiling
    //
    // Pair each sample with its matched control (if any), then gate on run_telomerehunter
    ch_th_input = ch_reads_validated
        .map { meta, _dt, bam, bai, _bed -> [ meta, bam, bai ] }
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

    // Append optional cytoband to fasta tuple for telomerehunter
    // Callers pass Channel.value([cytoband_path]) or Channel.value([[]])
    ch_fasta_cytoband = ch_fasta
        .combine(ch_cytoband)
        .map { meta, fasta, fai, cytoband -> [ meta, fasta, fai, cytoband ] }

    TELOMEREHUNTER(ch_th_input, ch_fasta_cytoband)

    ch_content_tsv          = TELOMEREHUNTER.out.summary
    ch_telomerehunter_tumor = TELOMEREHUNTER.out.tumor
    ch_telomerehunter_ctrl  = TELOMEREHUNTER.out.control

    //
    // Summarise: left-join length results with optional content results per sample
    //
    ch_summary_input = ch_length_with_tool
        .join(ch_content_tsv, remainder: true)
        .map { items ->
            def meta        = items[0]
            def length_tsv  = items[1]
            def tool        = items[2]
            def content_tsv = items.size() > 3 ? items[3] : []
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
