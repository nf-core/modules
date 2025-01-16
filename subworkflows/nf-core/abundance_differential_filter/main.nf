//
// Perform differential analysis
//

include { LIMMA_DIFFERENTIAL                  } from '../../../modules/nf-core/limma/differential/main'
include { LIMMA_DIFFERENTIAL as LIMMA_NORM    } from '../../../modules/nf-core/limma/differential/main'
include { DESEQ2_DIFFERENTIAL                 } from '../../../modules/nf-core/deseq2/differential/main'
include { DESEQ2_DIFFERENTIAL as DESEQ2_NORM  } from '../../../modules/nf-core/deseq2/differential/main'
include { PROPR_PROPD                         } from '../../../modules/nf-core/propr/propd/main'
include { CUSTOM_FILTERDIFFERENTIALTABLE      } from '../../../modules/nf-core/custom/filterdifferentialtable/main'

// Combine meta maps, including merging non-identical values of shared keys (e.g. 'id')
def mergeMaps(meta, meta2){
    (meta + meta2).collectEntries { k, v ->
        meta[k] && meta[k] != v ? [k, "${meta[k]}_${v}"] : [k, v]
    }
}

workflow ABUNDANCE_DIFFERENTIAL_FILTER {
    take:
    // Things we may need to iterate
    ch_input                 // [[meta_input], counts, analysis method, fc_threshold, stat_threshold]

    // Workflow-wide things, we don't need to iterate
    ch_samplesheet           // [ meta_exp, samplesheet ]
    ch_transcript_lengths    // [ meta_exp, transcript_lengths]
    ch_control_features      // [meta_exp, control_features]
    ch_contrasts             // [[ meta_contrast, contrast_variable, reference, target ]]

    main:

    // Set up how the channels crossed below will be used to generate channels for processing
    def criteria = multiMapCriteria { meta_input, abundance, analysis_method, fc_threshold, stat_threshold, meta_exp, samplesheet, meta_contrasts, variable, reference, target ->
        samples_and_matrix:
            meta_map = meta_input + [ 'method': analysis_method ]
            [meta_map, samplesheet, abundance]
        contrasts_for_diff:
            meta_map = mergeMaps(meta_contrasts, meta_input) + [ 'method': analysis_method ]
            [ meta_map, variable, reference, target ]
        filter_params:
            meta_map = mergeMaps(meta_contrasts, meta_input) + [ 'method': analysis_method ]
            [meta_map, [ 'fc_threshold': fc_threshold, 'stat_threshold': stat_threshold ]]
        contrasts_for_norm:
            meta_map = meta_input + [ 'method': analysis_method ]
            [ meta_map, variable, reference, target ]
    }

    // For DIFFERENTIAL modules we need to cross the things we're iterating so we
    // run differential analysis for every combination of matrix and contrast
    inputs = ch_input
        .combine(ch_samplesheet)
        .combine(ch_contrasts)
        .multiMap(criteria)

    // We only need a normalised matrix from one contrast. The reason we don't
    // simply use the first output from DIFFERENTIAL modules is that depending
    // on the contrast setting etc, these modules may subset matrices, hence
    // not returning the full normalized matrix as NORM modules would do.
    norm_inputs = ch_input
        .combine(ch_samplesheet)
        .combine(ch_contrasts.first()) // Just taking the first contrast
        .multiMap(criteria)

    inputs.samples_and_matrix.view{ it -> "inputs.samples_and_matrix is ${it}" }
    inputs.contrasts.view{ it -> "inputs.contrasts is ${it}" }
    inputs.filter_params.view{ it -> "inputs.filter_params is ${it}" }

    norm_inputs.samples_and_matrix.view{ it -> "norm_inputs.samples_and_matrix is ${it}" }
    norm_inputs.contrasts.view { it -> "norm_inputs.contrasts is ${it}" }
    norm_inputs.filter_params.view{ it -> "norm_inputs.filter_params is ${it}" }

    // ----------------------------------------------------
    // Run Limma
    // ----------------------------------------------------

    // NOTE that we run LIMMA_NORM just once to generate a normalised matrix.
    // As explained above, this is done to avoid obtaining a subset matrix
    // from LIMMA_DIFFERENTIAL.

    // Also NOTE that LIMMA_DIFFERENTIAL don't use the normalized matrix from
    // LIMMA_NORM directly. It internally runs normalization + DE analysis.

    LIMMA_NORM(
        norm_inputs.contrasts_for_norm.filter{it[0].method == 'limma'},
        norm_inputs.samples_and_matrix.filter{it[0].method == 'limma'}
    )

    LIMMA_DIFFERENTIAL(
        inputs.contrasts_for_diff.filter{ it[0].method == 'limma' },
        inputs.samples_and_matrix.filter{ it[0].method == 'limma' }
    )

    // ----------------------------------------------------
    // Run DESeq2
    // ----------------------------------------------------

    // NOTE that we run DESEQ2_NORM just once to generate a normalised matrix.
    // As explained above, this is done to avoid obtaining a subset matrix
    // from DESEQ2_DIFFERENTIAL.

    // Also NOTE that DESEQ2_DIFFERENTIAL don't use the normalized matrix from
    // DESEQ2_NORM directly. It internally runs normalization + DE analysis.

    DESEQ2_NORM(
        norm_inputs.contrasts_for_norm.filter{it[0].method == 'deseq2'},
        norm_inputs.samples_and_matrix.filter{it[0].method == 'deseq2'},
        ch_control_features.first(),
        ch_transcript_lengths.first()
    )

    DESEQ2_DIFFERENTIAL(
        inputs.contrasts_for_diff.filter{it[0].method == 'deseq2'},
        inputs.samples_and_matrix.filter{it[0].method == 'deseq2'},
        ch_control_features.first(),
        ch_transcript_lengths.first()
    )

    // ----------------------------------------------------
    // Run propd
    // ----------------------------------------------------

    // NOTE that this method don't rely on normalization, hence it does
    // not produce a normalized matrix.

    PROPR_PROPD(
        inputs.contrasts_for_diff.filter{it[0].method == 'propd'},
        inputs.samples_and_matrix.filter { it[0].method == 'propd' }
    )

    // ----------------------------------------------------
    // Collect results
    // ----------------------------------------------------

    ch_results = DESEQ2_DIFFERENTIAL.out.results
        .mix(LIMMA_DIFFERENTIAL.out.results)
        .mix(PROPR_PROPD.out.results_genewise)

    ch_normalised_matrix = DESEQ2_NORM.out.normalised_counts
        .mix(LIMMA_NORM.out.normalised_counts)

    ch_model = DESEQ2_DIFFERENTIAL.out.model
        .mix(LIMMA_DIFFERENTIAL.out.model)

    ch_versions = DESEQ2_DIFFERENTIAL.out.versions
        .mix(LIMMA_DIFFERENTIAL.out.versions)
        .mix(PROPR_PROPD.out.versions)

    // ----------------------------------------------------
    // Filter DE results
    // ----------------------------------------------------

    ch_diff_filter_params = ch_results
        .join(inputs.filter_params)
        .multiMap { meta, results, filter_meta ->
            def method_params = [
                'deseq2': [
                    fc_column: 'log2FoldChange', fc_cardinality: '>=',
                    stat_column: 'padj', stat_cardinality: '<='
                ],
                'limma' : [
                    fc_column: 'logFC', fc_cardinality: '>=',
                    stat_column: 'adj.P.Val', stat_cardinality: '<='
                ],
                'propd' : [
                    fc_column: 'lfc', fc_cardinality: '>=',
                    stat_column: 'weighted_connectivity', stat_cardinality: '>='
                ]
            ]
            filter_input: [meta + filter_meta, results]
            fc_input: [
                method_params[meta.method].fc_column,
                filter_meta.fc_threshold,
                method_params[meta.method].fc_cardinality
            ]
            stat_input: [
                method_params[meta.method].stat_column,
                filter_meta.stat_threshold,
                method_params[meta.method].stat_cardinality
            ]
        }

    // Filter differential results
    CUSTOM_FILTERDIFFERENTIALTABLE(
        ch_diff_filter_params.filter_input,
        ch_diff_filter_params.fc_input,
        ch_diff_filter_params.stat_input
    )

    ch_results.view{ it -> "ch_results is ${it}" }
    ch_normalised_matrix.view{ it -> "ch_normalised_matrix is ${it}" }
    ch_model.view{ it -> "ch_model is ${it}" }
    CUSTOM_FILTERDIFFERENTIALTABLE.out.filtered.view{ it -> "filtered output is ${it}" }

    emit:
    // main results
    results_genewise           = ch_results
    results_genewise_filtered  = CUSTOM_FILTERDIFFERENTIALTABLE.out.filtered

    // pairwise results
    adjacency                  = PROPR_PROPD.out.adjacency

    // other
    normalised_matrix          = ch_normalised_matrix
    variance_stabilised_matrix = DESEQ2_NORM.out.rlog_counts.mix(DESEQ2_NORM.out.vst_counts)
    model                      = ch_model
    versions                   = ch_versions
}
