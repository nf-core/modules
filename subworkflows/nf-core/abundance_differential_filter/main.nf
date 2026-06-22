//
// Perform differential analysis
//

include { LIMMA_DIFFERENTIAL                    } from '../../../modules/nf-core/limma/differential/main'
include { LIMMA_DIFFERENTIAL as LIMMA_NORM      } from '../../../modules/nf-core/limma/differential/main'
include { DESEQ2_DIFFERENTIAL                   } from '../../../modules/nf-core/deseq2/differential/main'
include { DESEQ2_DIFFERENTIAL as DESEQ2_NORM    } from '../../../modules/nf-core/deseq2/differential/main'
include { PROPR_PROPD                           } from '../../../modules/nf-core/propr/propd/main'
include { CUSTOM_FILTERDIFFERENTIALTABLE        } from '../../../modules/nf-core/custom/filterdifferentialtable/main'
include { VARIANCEPARTITION_DREAM               } from '../../../modules/nf-core/variancepartition/dream/main'
include { VARIANCEPARTITION_DREAM as DREAM_NORM } from '../../../modules/nf-core/variancepartition/dream/main'

// Combine meta maps, including merging non-identical values of shared keys (e.g. 'id')
def mergeMaps(meta, meta2){
    (meta + meta2).collectEntries { k, v ->
        meta[k] && meta[k] != v ? [k, "${meta[k]}_${v}"] : [k, v]
    }
}

def getDifferentialMethodParams(differential_method) {
    def method_params = [
        'deseq2': [
            differential_fc_column         : 'log2FoldChange',
            fc_cardinality                 : '>=',
            stat_column                    : 'padj',
            stat_cardinality               : '<='
        ],
        'limma' : [
            differential_fc_column         : 'logFC',
            fc_cardinality                 : '>=',
            stat_column                    : 'adj.P.Val',
            stat_cardinality               : '<='
        ],
        'propd' : [
            differential_fc_column         : 'LFC',
            fc_cardinality                 : '>=',
            stat_column                    : 'significant',
            stat_cardinality               : '<='
        ],
        'dream' : [
            differential_fc_column         : 'logFC',
            fc_cardinality                 : '>=',
            stat_column                    : 'adj.P.Val',
            stat_cardinality               : '<='
        ]
    ][differential_method]

    if (!method_params) {
        throw new IllegalArgumentException("Unsupported differential method: ${differential_method}")
    }

    method_params
}

workflow ABUNDANCE_DIFFERENTIAL_FILTER {
    take:
    // Things we may need to iterate
    ch_input                 // [ meta, abundance, analysis method, fc_threshold, stat_threshold ]

    // Workflow-wide things, we don't need to iterate
    ch_samplesheet           // [ meta, samplesheet ]
    ch_transcript_lengths    // [ meta, transcript_lengths]
    ch_control_features      // [ meta, control_features ]
    ch_contrasts             // [ meta, [meta_contrast], [variable], [reference], [target], [formula], [comparison] ]

    main:

    ch_versions = channel.empty()

    // Set up how the channels crossed below will be used to generate channels for processing
    def criteria = multiMapCriteria { meta, abundance, analysis_method, fc_threshold, stat_threshold, samplesheet, transcript_length, control_features, meta_contrast, variable, reference, target, formula, comparison ->
        def meta_with_method = meta + [ 'differential_method': analysis_method ]
        def meta_for_diff = mergeMaps(meta_contrast, meta_with_method)
        samples_and_matrix:
            [ meta_with_method, samplesheet, abundance ]
        contrasts_for_diff:
            [ meta_for_diff, variable, reference, target ]
        contrasts_for_diff_with_formula:
            [ meta_for_diff, variable, reference, target, formula, comparison ]
        filter_params:
            [ meta_for_diff, [ 'fc_threshold': fc_threshold, 'stat_threshold': stat_threshold ]]
        contrasts_for_norm:
            [ meta_with_method, variable, reference, target ]
        contrasts_for_norm_with_formula:
            [ meta_with_method, variable, reference, target, formula, comparison ]
        // these are optional files
        // return empty file if not available
        transcript_length:
            [ meta_with_method, (transcript_length)?:[] ]
        control_features:
            [ meta_with_method, (control_features)?:[] ]
    }

    // create joined channel with samplesheet and optional transcript lengths and control features
    ch_samplesheet_with_control = ch_samplesheet
        .join(ch_transcript_lengths, remainder:true)
        .join(ch_control_features, remainder:true)

    // For DIFFERENTIAL modules we need to cross the things we're iterating so we
    // run differential analysis for every combination of matrix and contrast
    inputs = ch_input
        .combine(ch_samplesheet_with_control, by:0)
        .combine(ch_contrasts.transpose(), by:0)
        .multiMap(criteria)

    // We only need a normalised matrix from one contrast. The reason we don't
    // simply use the first output from DIFFERENTIAL modules is that depending
    // on the contrast setting etc, these modules may subset matrices, hence
    // not returning the full normalized matrix as NORM modules would do.
    ch_first_contrast = ch_contrasts
        .map { meta, contrasts, variables, references, targets, formulas, comparisons ->
            tuple(meta, contrasts[0], variables[0], references[0], targets[0], formulas[0], comparisons[0])
        }
    norm_inputs = ch_input
        .combine(ch_samplesheet_with_control, by:0)
        .combine(ch_first_contrast, by:0)
        .multiMap(criteria)

    // ----------------------------------------------------
    // Run Limma
    // ----------------------------------------------------

    // NOTE that we run LIMMA_NORM just once to generate a normalised matrix.
    // As explained above, this is done to avoid obtaining a subset matrix
    // from LIMMA_DIFFERENTIAL.

    // Also NOTE that LIMMA_DIFFERENTIAL don't use the normalized matrix from
    // LIMMA_NORM directly. It internally runs normalization + DE analysis.

    LIMMA_NORM(
        norm_inputs.contrasts_for_norm_with_formula.filter{index -> index[0].differential_method == 'limma'},
        norm_inputs.samples_and_matrix.filter{index -> index[0].differential_method == 'limma'}
    )

    ch_versions = ch_versions.mix(LIMMA_NORM.out.versions.first())

    LIMMA_DIFFERENTIAL(
        inputs.contrasts_for_diff_with_formula.filter{index -> index[0].differential_method == 'limma' },
        inputs.samples_and_matrix.filter{index -> index[0].differential_method == 'limma' }
    )

    ch_versions = ch_versions.mix(LIMMA_DIFFERENTIAL.out.versions.first())

    // ----------------------------------------------------
    // Run DESeq2
    // ----------------------------------------------------

    // NOTE that we run DESEQ2_NORM just once to generate a normalised matrix.
    // As explained above, this is done to avoid obtaining a subset matrix
    // from DESEQ2_DIFFERENTIAL.

    // Also NOTE that DESEQ2_DIFFERENTIAL don't use the normalized matrix from
    // DESEQ2_NORM directly. It internally runs normalization + DE analysis.

    DESEQ2_NORM(
        norm_inputs.contrasts_for_norm_with_formula.filter{index -> index[0].differential_method == 'deseq2'},
        norm_inputs.samples_and_matrix.filter{index -> index[0].differential_method == 'deseq2'},
        norm_inputs.control_features.filter{index -> index[0].differential_method == 'deseq2'},
        norm_inputs.transcript_length.filter{index -> index[0].differential_method == 'deseq2'}
    )

    ch_versions = ch_versions.mix(DESEQ2_NORM.out.versions.first())

    DESEQ2_DIFFERENTIAL(
        inputs.contrasts_for_diff_with_formula.filter{index -> index[0].differential_method == 'deseq2'},
        inputs.samples_and_matrix.filter{index -> index[0].differential_method == 'deseq2'},
        inputs.control_features.filter{index -> index[0].differential_method == 'deseq2'},
        inputs.transcript_length.filter{index -> index[0].differential_method == 'deseq2'}
    )

    ch_versions = ch_versions.mix(DESEQ2_DIFFERENTIAL.out.versions.first())

    // ----------------------------------------------------
    // Run propd
    // ----------------------------------------------------

    // NOTE that this method don't rely on normalization, hence it does
    // not produce a normalized matrix.

    PROPR_PROPD(
        inputs.contrasts_for_diff.filter{index -> index[0].differential_method == 'propd'},
        inputs.samples_and_matrix.filter {index -> index[0].differential_method == 'propd' }
    )

    ch_versions = ch_versions.mix(PROPR_PROPD.out.versions.first())

    // ----------------------------------------------------
    // Run DREAM
    // ----------------------------------------------------

    // NOTE that we run DREAM_NORM just once to generate a normalised matrix.
    // As explained above, this is done to avoid obtaining a subset matrix
    // from VARIANCEPARTITION_DREAM.

    // Also NOTE that VARIANCEPARTITION_DREAM don't use the normalized matrix from
    // DREAM_NORM directly. It internally runs normalization + DE analysis.
    // Prepare DREAM inputs
    dream_inputs = inputs.contrasts_for_diff_with_formula
        .filter { meta, _variable, _reference, _target, _formula, _comparison ->
            meta.differential_method == 'dream'
        }

    DREAM_NORM(
        norm_inputs.contrasts_for_norm_with_formula.filter{index -> index[0].differential_method == 'dream'},
        norm_inputs.samples_and_matrix.filter{index -> index[0].differential_method == 'dream'}
    )

    VARIANCEPARTITION_DREAM(
        dream_inputs,
        inputs.samples_and_matrix.filter{index -> index[0].differential_method == 'dream' }
    )

    ch_versions = ch_versions.mix( VARIANCEPARTITION_DREAM.out.versions.first() )

    // ----------------------------------------------------
    // Collect results
    // ----------------------------------------------------

    ch_results = DESEQ2_DIFFERENTIAL.out.results
        .mix(LIMMA_DIFFERENTIAL.out.results)
        .mix(PROPR_PROPD.out.results_genewise)
        .mix(VARIANCEPARTITION_DREAM.out.results)

    ch_normalised_matrix = DESEQ2_NORM.out.normalised_counts
        .mix(LIMMA_NORM.out.normalised_counts)
        .mix(DREAM_NORM.out.normalised_counts)

    ch_model = DESEQ2_DIFFERENTIAL.out.model
        .mix(LIMMA_DIFFERENTIAL.out.model)
        .mix(VARIANCEPARTITION_DREAM.out.model)

    ch_variance_stabilised_matrix = DESEQ2_NORM.out.rlog_counts
        .mix(DESEQ2_NORM.out.vst_counts)
        .groupTuple()

    ch_dispersion_plot = DESEQ2_DIFFERENTIAL.out.dispersion_plot_png

    ch_md_plot = LIMMA_DIFFERENTIAL.out.md_plot

    ch_rdata = DESEQ2_DIFFERENTIAL.out.rdata
        .mix(LIMMA_DIFFERENTIAL.out.rdata)

    ch_size_factors = DESEQ2_DIFFERENTIAL.out.size_factors

    ch_session_info = DESEQ2_DIFFERENTIAL.out.session_info
        .mix(LIMMA_DIFFERENTIAL.out.session_info)

    // ----------------------------------------------------
    // Filter DE results
    // ----------------------------------------------------

    ch_diff_filter_params = ch_results
        .join(inputs.filter_params)
        .multiMap { meta, results, filter_meta ->
            def method_params = getDifferentialMethodParams(meta.differential_method)
            filter_input: [meta + filter_meta, results]
            fc_input: [
                method_params.differential_fc_column,
                filter_meta.fc_threshold,
                method_params.fc_cardinality
            ]
            stat_input: [
                method_params.stat_column,
                filter_meta.stat_threshold,
                method_params.stat_cardinality
            ]
        }

    // Filter differential results
    CUSTOM_FILTERDIFFERENTIALTABLE(
        ch_diff_filter_params.filter_input,
        ch_diff_filter_params.fc_input,
        ch_diff_filter_params.stat_input
    )

    emit:

    // main results
    results_genewise           = ch_results
    results_genewise_filtered  = CUSTOM_FILTERDIFFERENTIALTABLE.out.filtered

    // pairwise results
    adjacency                  = PROPR_PROPD.out.adjacency

    // other matrices
    normalised_matrix          = ch_normalised_matrix
    variance_stabilised_matrix = ch_variance_stabilised_matrix

    // differential model
    model                      = ch_model

    // size factors
    size_factors               = ch_size_factors

    // plots
    dispersion_plot            = ch_dispersion_plot
    md_plot                    = ch_md_plot

    // rdata and session info
    rdata                      = ch_rdata
    session_info               = ch_session_info

    // versions
    versions                   = ch_versions
}
