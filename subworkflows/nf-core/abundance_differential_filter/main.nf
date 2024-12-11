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
    ch_abundance          // [ meta_exp, counts ] with meta keys: method, args_diff
    ch_transcript_lengths // [ meta_exp, transcript_lengths]
    ch_control_features   // [ meta_exp, control_features]
    ch_samplesheet        // [ meta_exp, samplesheet ]
    ch_contrasts          // [ meta_contrast, contrast_variable, reference, target ]
    differential_method   // limma, deseq2, propd
    FC_threshold          // float
    padj_threshold        // float

    main:

    // initialize empty results channels
    ch_normalised_matrix           = Channel.empty()
    ch_variance_stabilised_matrix  = Channel.empty()
    ch_model                       = Channel.empty()

    ch_results_genewise            = Channel.empty()
    ch_results_genewise_filtered   = Channel.empty()
    ch_versions                    = Channel.empty()

    // Derive some commonly used derived channels
    ch_samples_and_matrix = ch_samplesheet.join(ch_abundance).first()

    logFC_column = null
    padj_column = null

    // ----------------------------------------------------
    // Perform differential analysis with DESeq2
    // ----------------------------------------------------

    if (differential_method == 'deseq2') {

        logFC_column = "log2FoldChange"
        padj_column  = "padj"

        // Run the module once to generate a normalised matrix. We can't just
        // use e.g. the first run of DESEQ_DIFFERENTIAL, because it may remove
        // some samples

        DESEQ2_NORM (
            ch_contrasts.first(),
            ch_samples_and_matrix,
            ch_control_features,
            ch_transcript_lengths
        )

        ch_normalised_matrix = DESEQ2_NORM.out.normalised_counts
        ch_variance_stabilised_matrix = DESEQ2_NORM.out.rlog_counts.concat(DESEQ2_NORM.out.vst_counts)

        // Run the DESeq differential module

        DESEQ2_DIFFERENTIAL (
            ch_contrasts,
            ch_samples_and_matrix,
            ch_control_features,
            ch_transcript_lengths
        )

        ch_results_genewise = DESEQ2_DIFFERENTIAL.out.results
        ch_model = DESEQ2_DIFFERENTIAL.out.model

        ch_versions = ch_versions
            .mix(DESEQ2_DIFFERENTIAL.out.versions)

    } else if (differential_method == 'limma'){

        logFC_column = "logFC"
        padj_column  = "adj.P.Val"

        // Run the module once to generate a normalised matrix. We can't just
        // use e.g. the first run of LIMMA_DIFFERENTIAL, because it may remove
        // some samples

        LIMMA_NORM (
            ch_contrasts.first(),
            ch_samples_and_matrix
        )

        ch_normalised_matrix = LIMMA_NORM.out.normalised_counts

        LIMMA_DIFFERENTIAL (
            ch_contrasts,
            ch_samples_and_matrix
        )
        ch_results_genewise = LIMMA_DIFFERENTIAL.out.results
        ch_model = LIMMA_DIFFERENTIAL.out.model

        ch_versions = ch_versions
            .mix(LIMMA_DIFFERENTIAL.out.versions)

    }

    CUSTOM_FILTERDIFFERENTIALTABLE(
        ch_results_genewise,
        logFC_column,
        FC_threshold,
        padj_column,
        padj_threshold
    )

    emit:
    results_genewise           = ch_results_genewise               // channel: [ val(meta), path(results) ]
    results_genewise_filtered  = CUSTOM_FILTERDIFFERENTIALTABLE.out.filtered // channel: [ val(meta), path(filtered_results) ]
    normalised_matrix          = ch_normalised_matrix              // channel: [ val(meta), path(normalised_matrix) ]
    variance_stabilised_matrix = ch_variance_stabilised_matrix     // channel: [ val(meta), path(variance_stabilised_matrix) ] (Optional)
    model                      = ch_model                          // channel: [ val(meta), path(model) ]
    versions                   = ch_versions                       // channel: [ path(versions.yml) ]
}
