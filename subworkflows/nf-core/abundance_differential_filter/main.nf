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
    ch_input                 // [[meta_input], counts, analysis method,]

    // Workflow-wide things, we don't need to iterate
    ch_samplesheet           // [ meta_exp, samplesheet ]
    ch_transcript_lengths    // [ meta_exp, transcript_lengths]
    ch_control_features      // [meta_exp, control_features]
    ch_contrasts             // [[ meta_contrast, contrast_variable, reference, target ]]
    FC_threshold             // FC threshold (value)
    padj_threshold           // padj threshold (value)


    main:

    // We need to cross the things we're iterating
    inputs = ch_input
        .combine(ch_samplesheet)
        .combine(ch_contrasts)
        .multiMap { meta_input, abundance, analysis_method, meta_exp, samplesheet, contrast_meta, variable, reference, target ->
            samples_and_matrix: [meta_input + [ 'method': analysis_method ], samplesheet, abundance]
            contrasts: [ contrast_meta + [ 'method': analysis_method ], variable, reference, target]
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

    // Define logFC and padj columns based on method in ch_results
    ch_columns = ch_results.map { meta, results ->
        def columns = meta.method == 'deseq2' ? ["log2FoldChange", "padj"] : ["logFC", "adj.P.Val"]
        [meta, columns]
    }

    // Filter differential results
    CUSTOM_FILTERDIFFERENTIALTABLE(
        ch_results,
        ch_columns.map { meta, columns -> columns[0] },  // logFC column
        FC_threshold,
        ch_columns.map { meta, columns -> columns[1] },  // padj column
        padj_threshold
    )

    emit:
    results_genewise           = ch_results
    results_genewise_filtered  = CUSTOM_FILTERDIFFERENTIALTABLE.out.filtered
    normalised_matrix          = ch_normalised_matrix
    variance_stabilised_matrix = DESEQ2_NORM.out.rlog_counts.mix(DESEQ2_NORM.out.vst_counts)
    model                      = ch_model
    versions                   = ch_versions
}
