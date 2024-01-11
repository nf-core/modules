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
    ch_versions = Channel.empty()

    // MODULE: SEQKIT_SEQ as FILTER_BY_LENGTH
    FILTER_BY_LENGTH ( ch_fasta )

    ch_filtered_fasta       = FILTER_BY_LENGTH.out.fastx
    ch_versions             = ch_versions.mix(FILTER_BY_LENGTH.out.versions.first())

    // MODULE: SEQKIT_SORT as SORT_BY_LENGTH
    SORT_BY_LENGTH ( ch_filtered_fasta )

    ch_sorted_fasta         = SORT_BY_LENGTH.out.fastx
    ch_versions             = ch_versions.mix(SORT_BY_LENGTH.out.versions.first())

    // TIDK_EXPLORE
    TIDK_EXPLORE ( ch_filtered_fasta )

    ch_top_sequence         = TIDK_EXPLORE.out.top_sequence
    ch_versions             = ch_versions.mix(TIDK_EXPLORE.out.versions.first())

    // TIDK_SEARCH as TIDK_SEARCH_APRIORI
    ch_apriori_inputs       = ch_sorted_fasta
                            | map { meta, fasta -> [ meta.id, meta, fasta ] }
                            | join(
                                ( ch_apriori_sequence ?: Channel.empty() )
                                | map { meta, seq -> [ meta.id, seq ] }
                            )
                            | map { id, meta, fasta, seq -> [ meta, fasta, seq ] }

    TIDK_SEARCH_APRIORI (
        ch_apriori_inputs.map { meta, fasta, seq -> [ meta, fasta ] },
        ch_apriori_inputs.map { meta, fasta, seq -> seq }
    )

    ch_apriori_tsv          = TIDK_SEARCH_APRIORI.out.tsv
    ch_versions             = ch_versions.mix(TIDK_SEARCH_APRIORI.out.versions.first())

    // TIDK_SEARCH as TIDK_SEARCH_APOSTERIORI
    ch_aposteriori_inputs   = ch_sorted_fasta
                            | join(ch_top_sequence)
                            | map { meta, fasta, txt ->
                                [ meta, fasta, txt.getText().strip() ]
                            }

    TIDK_SEARCH_APOSTERIORI (
        ch_aposteriori_inputs.map { meta, fasta, seq -> [ meta, fasta ] },
        ch_aposteriori_inputs.map { meta, fasta, seq -> seq }
    )

    ch_aposteriori_tsv      = TIDK_SEARCH_APOSTERIORI.out.tsv
    ch_versions             = ch_versions.mix(TIDK_SEARCH_APOSTERIORI.out.versions.first())

    // TIDK_PLOT as TIDK_PLOT_APRIORI
    TIDK_PLOT_APRIORI ( ch_apriori_tsv )

    ch_apriori_svg          = TIDK_PLOT_APRIORI.out.svg
    ch_versions             = ch_versions.mix(TIDK_PLOT_APRIORI.out.versions.first())

    // TIDK_PLOT as TIDK_PLOT_APOSTERIORI
    TIDK_PLOT_APOSTERIORI ( ch_aposteriori_tsv )

    ch_aposteriori_svg      = TIDK_PLOT_APOSTERIORI.out.svg
    ch_versions             = ch_versions.mix(TIDK_PLOT_APOSTERIORI.out.versions.first())

    emit:
    apriori_tsv             = ch_apriori_tsv        // channel: [ val(meta), tsv ]
    apriori_svg             = ch_apriori_svg        // channel: [ val(meta), svg ]
    aposteriori_sequence    = ch_top_sequence       // channel: [ val(meta), txt ]
    aposteriori_tsv         = ch_aposteriori_tsv    // channel: [ val(meta), tsv ]
    aposteriori_svg         = ch_aposteriori_svg    // channel: [ val(meta), svg ]
    versions                = ch_versions           // channel: [ versions.yml ]
}
