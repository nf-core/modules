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
    ch_abundance             // [[ meta_exp, counts ]] with meta keys: method, args_diff
    ch_contrasts             // [[ meta_contrast, contrast_variable, reference, target ]]
    ch_analysis_method       // [[ 'limma|deseq2' ]]
    ch_samplesheet           // [ meta_exp, samplesheet ]
    ch_transcript_lengths    // [ meta_exp, transcript_lengths]
    ch_control_features      // [meta_exp, control_features]
    FC_threshold             // float
    padj_threshold           // float

    main:

    // Combine and prepare data for differential analysis
    // This operation creates all possible combinations of abundance data, contrasts, and analysis methods
    // It also incorporates sample information and adds method details to each contrast

    ch_combos = ch_samplesheet
        .combine(ch_abundance)
        .map{ samples_meta, samples, abundance_meta, abundance -> [ samples_meta + abundance_meta, samples, abundance ] }
        .merge(ch_contrasts) { abundance, contrasts -> [abundance, contrasts] }
        .merge(ch_analysis_method) { prev, methods -> [ prev[0], prev[1], methods ] }
        .flatMap { samples_abundance, contrasts, methods ->
            contrasts.collectMany { contrast ->
                methods.collect { method ->
                    [
                        samples_abundance: samples_abundance,
                        contrast: [
                            contrast[0] + [method: method]  // Add method to the map in the first element
                        ] + contrast.drop(1)  // Add all remaining elements of the contrast tuple
                    ]
                }
            }
        }

    ch_samples_and_matrix = ch_combos.map{it.samples_abundance}
    ch_contrasts = ch_combos.map{it.contrast}

    // Perform normalization and differential analysis
    DESEQ2_NORM(
        ch_contrasts.filter { it[0].method == 'deseq2' }.first(),
        ch_samples_and_matrix,
        ch_control_features,
        ch_transcript_lengths
    )

    LIMMA_NORM(
        ch_contrasts.filter { it[0].method == 'limma' }.first(),
        ch_samples_and_matrix
    )

    DESEQ2_DIFFERENTIAL(
        ch_contrasts.filter { it[0].method == 'deseq2' }.first(),
        ch_samples_and_matrix,
        ch_control_features,
        ch_transcript_lengths
    )

    LIMMA_DIFFERENTIAL(
        ch_contrasts.filter { it[0].method == 'limma' }.first(),
        ch_samples_and_matrix
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
        Channel.value(FC_threshold),
        ch_columns.map { meta, columns -> columns[1] },  // padj column
        Channel.value(padj_threshold)
    )

    emit:
    results_genewise           = ch_results
    results_genewise_filtered  = CUSTOM_FILTERDIFFERENTIALTABLE.out.filtered
    normalised_matrix          = ch_normalised_matrix
    variance_stabilised_matrix = DESEQ2_NORM.out.rlog_counts.mix(DESEQ2_NORM.out.vst_counts)
    model                      = ch_model
    versions                   = ch_versions
}
