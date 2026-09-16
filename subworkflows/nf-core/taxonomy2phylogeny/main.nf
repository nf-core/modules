include { CUSTOM_TAXONOMYTREE } from '../../../modules/nf-core/custom/taxonomytree/main'
include { RAXMLNG_SEARCH      } from '../../../modules/nf-core/raxmlng/search/main'

workflow TAXONOMY2PHYLOGENY {

    take:
    ch_taxonomy_alignment // channel: [ val(meta), path(taxonomy), path(alignment), val(raxmlng_model) ]

    main:
    CUSTOM_TAXONOMYTREE(
        ch_taxonomy_alignment.map { meta, taxonomy, alignment, raxmlng_model -> [ meta, taxonomy ] }
    )

    // RAXMLNG_SEARCH's tree_constraint is a bare path, not a [meta, path] tuple, so it
    // can't be paired with the right alignment by a plain .map() off CUSTOM_TAXONOMYTREE's
    // own output -- .join() first, on the shared meta, then re-project the single joined
    // channel into RAXMLNG_SEARCH's separate positional inputs, so both stay correctly
    // paired regardless of the two processes' relative completion order.
    def ch_search_input = ch_taxonomy_alignment
        .map { meta, taxonomy, alignment, raxmlng_model -> [ meta, alignment, raxmlng_model ] }
        .join(CUSTOM_TAXONOMYTREE.out.guide_tree)
    // ch_search_input: [ meta, alignment, raxmlng_model, guide_tree ]

    RAXMLNG_SEARCH(
        ch_search_input.map { meta, alignment, raxmlng_model, _guide_tree -> [ meta, alignment, raxmlng_model ] },
        [],
        ch_search_input.map { meta, alignment, raxmlng_model, guide_tree -> guide_tree },
        [],
        []
    )

    emit:
    tree       = RAXMLNG_SEARCH.out.phylogeny         // channel: [ val(meta), path(tree) ]      ML phylogeny, consistent with the taxonomy constraint
    model      = RAXMLNG_SEARCH.out.best_model        // channel: [ val(meta), path(model) ]     RAxML-NG model file, consistent with tree + alignment
    guide_tree = CUSTOM_TAXONOMYTREE.out.guide_tree   // channel: [ val(meta), path(guide_tree) ] the multifurcating taxonomy-only guide tree
}
