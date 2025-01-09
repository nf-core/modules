//
// Perform differential analysis
//

include { LIMMA_DIFFERENTIAL                  } from '../../../modules/nf-core/limma/differential/main'
include { LIMMA_DIFFERENTIAL as LIMMA_NORM    } from '../../../modules/nf-core/limma/differential/main'
include { DESEQ2_DIFFERENTIAL                 } from '../../../modules/nf-core/deseq2/differential/main'
include { DESEQ2_DIFFERENTIAL as DESEQ2_NORM  } from '../../../modules/nf-core/deseq2/differential/main'
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
    ch_input                 // [[meta_input], counts, analysis method, fc_threshold, padj_threshold]

    // Workflow-wide things, we don't need to iterate
    ch_samplesheet           // [ meta_exp, samplesheet ]
    ch_transcript_lengths    // [ meta_exp, transcript_lengths]
    ch_control_features      // [meta_exp, control_features]
    ch_contrasts             // [[ meta_contrast, contrast_variable, reference, target ]]

    main:

    // Set up how the channels crossed below will be used to generate channels for processing
    def criteria = multiMapCriteria { meta_input, abundance, analysis_method, fc_threshold, padj_threshold, meta_exp, samplesheet, meta_contrasts, variable, reference, target ->
        samples_and_matrix:
            meta_map = meta_input + [ 'method': analysis_method ]
            [meta_map, samplesheet, abundance]
        contrasts:
            meta_map = mergeMaps(meta_contrasts, meta_input) + [ 'method': analysis_method ]
            [ meta_map, variable, reference, target ]
        filter_params:
            meta_map = mergeMaps(meta_contrasts, meta_input) + [ 'method': analysis_method ]
            [meta_map, [ 'fc_threshold': fc_threshold, 'padj_threshold': padj_threshold ]]
    }

    // For differential we need to cross the things we're iterating so we run
    // differential analysis for every combination of matrix and contrast
    inputs = ch_input
        .combine(ch_samplesheet)
        .combine(ch_contrasts)
        .multiMap(criteria)

    // We only need a normalised matrix from one contrast. The reason we don't
    //just use the output from the first differential is that the methods can
    // subset matrices
    norm_inputs = ch_input
        .combine(ch_samplesheet)
        .combine(ch_contrasts.first()) // Just taking the first contrast
        .multiMap(criteria)

    // Perform normalization and differential analysis
    DESEQ2_NORM(
        norm_inputs.contrasts.filter{it[0].method == 'deseq2'}.first(),
        norm_inputs.samples_and_matrix.filter{it[0].method == 'deseq2'},
        ch_control_features.first(),
        ch_transcript_lengths.first()
    )

    LIMMA_NORM(
        norm_inputs.contrasts.filter{it[0].method == 'limma'}.first(),
        norm_inputs.samples_and_matrix.filter{it[0].method == 'limma'}
    )

    DESEQ2_DIFFERENTIAL(
        inputs.contrasts.filter{it[0].method == 'deseq2'},
        inputs.samples_and_matrix.filter{it[0].method == 'deseq2'},
        ch_control_features.first(),
        ch_transcript_lengths.first()
    )

    LIMMA_DIFFERENTIAL(
        inputs.contrasts.filter{it[0].method == 'limma'},
        inputs.samples_and_matrix.filter { it[0].method == 'limma' }
    )

    // Combine results
    ch_results = DESEQ2_DIFFERENTIAL.out.results.mix(LIMMA_DIFFERENTIAL.out.results)
    ch_normalised_matrix = DESEQ2_NORM.out.normalised_counts.mix(LIMMA_NORM.out.normalised_counts)
    ch_model = DESEQ2_DIFFERENTIAL.out.model.mix(LIMMA_DIFFERENTIAL.out.model)
    ch_versions = DESEQ2_DIFFERENTIAL.out.versions
        .mix(LIMMA_DIFFERENTIAL.out.versions)

    // Extract the fc and pval filters from the metamap we stashed them in
    ch_diff_filter_params = ch_results
        .join(inputs.filter_params)
        .multiMap { meta, results, filter_meta ->
            filter_input: [meta + filter_meta, results]
            fc_column: meta.method == 'deseq2' ? 'log2FoldChange' : 'logFC'
            padj_column: meta.method == 'deseq2' ? 'padj' : 'adj.P.Val'
            fc_threshold: filter_meta.fc_threshold
            padj_threshold: filter_meta.padj_threshold
        }

    // Filter differential results
    CUSTOM_FILTERDIFFERENTIALTABLE(
        ch_diff_filter_params.filter_input,
        ch_diff_filter_params.fc_column,
        ch_diff_filter_params.fc_threshold,
        ch_diff_filter_params.padj_column,
        ch_diff_filter_params.padj_threshold
    )

    emit:
    results_genewise           = ch_results.map{meta, results -> [meta, results] }
    results_genewise_filtered  = CUSTOM_FILTERDIFFERENTIALTABLE.out.filtered
    normalised_matrix          = ch_normalised_matrix
    variance_stabilised_matrix = DESEQ2_NORM.out.rlog_counts.mix(DESEQ2_NORM.out.vst_counts)
    model                      = ch_model
    versions                   = ch_versions
}
