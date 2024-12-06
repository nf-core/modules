//
// Perform differential analysis
//

include { LIMMA_DIFFERENTIAL                  } from '../../../modules/nf-core/limma/differential/main'
include { LIMMA_DIFFERENTIAL as LIMMA_NORM    } from '../../../modules/nf-core/limma/differential/main'
include { DESEQ2_DIFFERENTIAL                 } from '../../../modules/nf-core/deseq2/differential/main'
include { DESEQ2_DIFFERENTIAL as DESEQ2_NORM  } from '../../../modules/nf-core/deseq2/differential/main'
include { CUSTOM_FILTERDIFFERENTIALTABLE      } from '../../../modules/nf-core/custom/filterdifferentialtable/main'

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

    // We need to cross the things we're iterating
    inputs = ch_input
        .combine(ch_samplesheet)
        .combine(ch_contrasts)
        .multiMap { meta_input, abundance, analysis_method, fc_threshold, padj_threshold, meta_exp, samplesheet, contrast_meta, variable, reference, target ->
            samples_and_matrix:
                meta_map = meta_input + [ 'method': analysis_method, 'fc': fc_threshold, 'padj': padj_threshold ]
                [meta_map, samplesheet, abundance]
            contrasts:
                meta_map = contrast_meta + [ 'method': analysis_method ]
                [ meta_map, variable, reference, target]

        }

    // Perform normalization and differential analysis
    DESEQ2_NORM(
        inputs.contrasts.filter{it[0].method == 'deseq2'}.first(),
        inputs.samples_and_matrix.filter{it[0].method == 'deseq2'},
        ch_control_features,
        ch_transcript_lengths
    )

    LIMMA_NORM(
        inputs.contrasts.filter{it[0].method == 'limma'}.first(),
        inputs.samples_and_matrix.filter{it[0].method == 'limma'}
    )

    DESEQ2_DIFFERENTIAL(
        inputs.contrasts.filter{it[0].method == 'deseq2'},
        inputs.samples_and_matrix.filter{it[0].method == 'deseq2'},
        ch_control_features,
        ch_transcript_lengths
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
    ch_diff_filter_params = ch_results.multiMap { meta, results ->
        fc_column: meta.method == 'deseq2' ? 'log2FoldChange' : 'logFC'
        padj_column: meta.method == 'deseq2' ? 'padj' : 'adj.P.Val'
        fc_threshold: meta.fc
        padj_threshold: meta.padj
    }

    // Filter differential results
    CUSTOM_FILTERDIFFERENTIALTABLE(
        ch_results,
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
