include { OPT_FLIP  } from '../../../modules/nf-core/opt/flip/main'
include { OPT_TRACK } from '../../../modules/nf-core/opt/track/main'
include { OPT_STAT  } from '../../../modules/nf-core/opt/stat/main'


workflow OPT_FLIP_TRACK_STAT {
    take:
    ch_probe_fasta   // channel: [ val(meta), [ "panel_probes_sequences.fasta" ] ]
    ch_references    // channel: [ val(meta), ["reference_annotations.gff"], ["reference_annotations.fa"] ]
    ch_gene_synonyms // channel: [ "path-to-gene-synonyms" ]

    main:

    ch_summary = channel.empty()

    // correct probes that are aligning to opposite strand with `flip`
    OPT_FLIP(ch_probe_fasta, ch_references)

    // align query probe sequences to target transcriptome
    OPT_TRACK(OPT_FLIP.out.fwd_oriented_fa, ch_references)

    // summarizes opt binding predictions
    OPT_STAT(OPT_TRACK.out.probes2target, OPT_FLIP.out.fwd_oriented_fa, ch_gene_synonyms)

    ch_summary = OPT_STAT.out.summary

    emit:
    summary  = ch_summary  // channel: [ val(meta), ["collapsed_summary.tsv"]]
}