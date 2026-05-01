include { SEQKIT_SEQ as FILTER_BY_LENGTH            } from '../../../modules/nf-core/seqkit/seq'
include { SEQKIT_SORT as SORT_BY_LENGTH             } from '../../../modules/nf-core/seqkit/sort'
include { TIDK_EXPLORE                              } from '../../../modules/nf-core/tidk/explore'
include { TIDK_SEARCH as TIDK_SEARCH_APRIORI        } from '../../../modules/nf-core/tidk/search'
include { TIDK_SEARCH as TIDK_SEARCH_APOSTERIORI    } from '../../../modules/nf-core/tidk/search'
include { TIDK_PLOT as TIDK_PLOT_APRIORI            } from '../../../modules/nf-core/tidk/plot'
include { TIDK_PLOT as TIDK_PLOT_APOSTERIORI        } from '../../../modules/nf-core/tidk/plot'


workflow FASTA_EXPLORE_SEARCH_PLOT_TIDK {

    take:
    ch_fasta                // channel: [ val(meta), [ fasta ] ]
    ch_apriori_sequence     // channel: [ val(meta), val(sequence) ]; Optional: Set to [] if not needed
                            // val(meta) from ch_fasta and ch_apriori_sequence are only required to have
                            // the same `id`

    main:

    // MODULE: SEQKIT_SEQ as FILTER_BY_LENGTH
    FILTER_BY_LENGTH ( ch_fasta )

    ch_filtered_fasta       = FILTER_BY_LENGTH.out.fastx

    // MODULE: SEQKIT_SORT as SORT_BY_LENGTH
    SORT_BY_LENGTH ( ch_filtered_fasta )

    ch_sorted_fasta         = SORT_BY_LENGTH.out.fastx

    // TIDK_EXPLORE
    TIDK_EXPLORE ( ch_filtered_fasta )

    ch_top_sequence         = TIDK_EXPLORE.out.top_sequence

    // TIDK_SEARCH as TIDK_SEARCH_APRIORI
    ch_apriori_inputs       = ch_sorted_fasta
                            | map { meta, fasta -> [ meta.id, meta, fasta ] }
                            | combine(
                                ( ch_apriori_sequence ?: channel.empty() )
                                | map { meta, seq -> [ meta.id, seq ] },
                                by:0
                            )
                            | map { _id, meta, fasta, seq -> [ meta, fasta, seq ] }

    ch_apriori_inputs.view()

    TIDK_SEARCH_APRIORI (
        ch_apriori_inputs.map { meta, fasta, _seq -> [ meta, fasta ] },
        ch_apriori_inputs.map { _meta, _fasta, seq -> seq }
    )

    ch_apriori_tsv          = TIDK_SEARCH_APRIORI.out.tsv

    // TIDK_SEARCH as TIDK_SEARCH_APOSTERIORI
    ch_aposteriori_inputs   = ch_sorted_fasta
                            | combine(ch_top_sequence, by:0)
                            | map { meta, fasta, txt ->
                                [ meta, fasta, txt.getText().strip() ]
                            }

    TIDK_SEARCH_APOSTERIORI (
        ch_aposteriori_inputs.map { meta, fasta, _seq -> [ meta, fasta ] },
        ch_aposteriori_inputs.map { _meta, _fasta, seq -> seq }
    )

    ch_aposteriori_tsv      = TIDK_SEARCH_APOSTERIORI.out.tsv

    // TIDK_PLOT as TIDK_PLOT_APRIORI
    TIDK_PLOT_APRIORI ( ch_apriori_tsv )

    ch_apriori_svg          = TIDK_PLOT_APRIORI.out.svg

    // TIDK_PLOT as TIDK_PLOT_APOSTERIORI
    TIDK_PLOT_APOSTERIORI ( ch_aposteriori_tsv )

    ch_aposteriori_svg      = TIDK_PLOT_APOSTERIORI.out.svg

    emit:
    apriori_tsv             = ch_apriori_tsv        // channel: [ val(meta), tsv ]
    apriori_svg             = ch_apriori_svg        // channel: [ val(meta), svg ]
    aposteriori_sequence    = ch_top_sequence       // channel: [ val(meta), txt ]
    aposteriori_tsv         = ch_aposteriori_tsv    // channel: [ val(meta), tsv ]
    aposteriori_svg         = ch_aposteriori_svg    // channel: [ val(meta), svg ]
}
