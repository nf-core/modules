
//
// Perform enrichment analysis
//
include { GPROFILER2_GOST          } from '../../../modules/nf-core/gprofiler2/gost/main.nf'
include { CUSTOM_TABULARTOGSEAGCT  } from '../../../modules/nf-core/custom/tabulartogseagct/main.nf'
include { CUSTOM_TABULARTOGSEACLS  } from '../../../modules/nf-core/custom/tabulartogseacls/main.nf'
include { CUSTOM_TABULARTOGSEACHIP } from '../../../modules/nf-core/custom/tabulartogseachip/main.nf'
include { GSEA_GSEA                } from '../../../modules/nf-core/gsea/gsea/main.nf'
include { PROPR_GREA               } from '../../../modules/nf-core/propr/grea/main.nf'
include { DECOUPLER_DECOUPLER      } from '../../../modules/nf-core/decoupler/decoupler/main'

// Combine meta maps, including merging non-identical values of shared keys (e.g. 'id')
def mergeMaps(meta, meta2){
    (meta + meta2).collectEntries { k, v ->
        meta[k] && meta[k] != v ? [k, "${meta[k]}_${v}"] : [k, v]
    }
}

workflow DIFFERENTIAL_FUNCTIONAL_ENRICHMENT {
    take:
    // input data for functional analysis
    // Note that genesets and background are optional depending on the method.
    // Please set to [] if not provided, eg: [meta, input, [], [], method]
    ch_input                 // [ meta, input file, genesets file, background file, method to run ]

    // other - for the moment these files are only needed for GSEA
    // as it is the only one that takes expression data as input
    // if in the future this setting is changed, this section could be removed
    ch_contrasts             // [ meta, [meta_contrast], [variable], [reference], [target], [formula], [comparison] ]
    ch_samplesheet           // [ meta, samples sheet ]
    ch_featuresheet          // [ meta, features sheet, features id, features symbol ]

    main:

    ch_versions = channel.empty()

    // Add method information into meta map of ch_input
    // This information is used later to determine which method to run for each input
    // Also, reorganize the structure to match them with the modules' input organization
    ch_input_for_other = ch_input
        .join(ch_featuresheet)
        .multiMap { meta, file, genesets, background, method, features_sheet, features_id, features_symbol ->
            def meta_with_method = meta + [ 'functional_method': method ]
            input:
                [ meta_with_method, file ]
            genesets:
                [ meta_with_method, genesets ]
            background:
                [ meta_with_method, background ]
            features:
                [ meta_with_method, features_sheet, features_id, features_symbol]
        }

    // In the case of GSEA, it needs additional files coming from other channels that other methods don't use
    // here we define the input channel for the GSEA section

    def criteria = multiMapCriteria { meta, input, genesets, _background, method, samplesheet, featuresheet, features_id, features_symbol, meta_contrast, _variable, _reference, _target, _formula, _comparison ->
        def meta_with_method = meta + [ 'functional_method': method ]
        def meta_with_contrast = mergeMaps(meta_contrast, meta_with_method)
        input:
            [ meta_with_contrast, input ]
        genesets:
            [ meta_with_contrast, genesets ]
        contrasts_and_samples:
            [ meta_with_contrast, samplesheet ]
        features:
            [ meta_with_method, featuresheet ]
        features_cols:
            [ features_id, features_symbol ]
    }

    // GSEA uses meta.variable, so only keep contrasts where meta.variable is present
    ch_contrasts_transposed = ch_contrasts.transpose()
        .filter { _meta, _contrastMap, variable, _reference, _target, _formula, _comparison ->
            variable?.trim()
        }

    ch_input_for_gsea = ch_input
        .filter{ it -> it[4] == 'gsea' }
        .combine(ch_samplesheet.join(ch_featuresheet), by:0)
        .combine(ch_contrasts_transposed, by:0)
        .multiMap(criteria)


    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    GPROFILER2_GOST(
        ch_input_for_other.input.filter{ index -> index[0].functional_method == 'gprofiler2' },
        ch_input_for_other.genesets.filter{ index -> index[0].functional_method == 'gprofiler2'},
        ch_input_for_other.background.filter{ index -> index[0].functional_method == 'gprofiler2'}
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // NOTE that GCT input can be more than 1, if they come from different tools (eg. limma, deseq2).
    // CLS input can be as many as combinations of input x contrasts
    // Whereas features can be only one file.

    CUSTOM_TABULARTOGSEAGCT(ch_input_for_gsea.input)

    CUSTOM_TABULARTOGSEACLS(ch_input_for_gsea.contrasts_and_samples)

    CUSTOM_TABULARTOGSEACHIP(
        ch_input_for_gsea.features.first(),
        ch_input_for_gsea.features_cols.first()
    )

    ch_input_for_gsea = CUSTOM_TABULARTOGSEAGCT.out.gct
        .join(CUSTOM_TABULARTOGSEACLS.out.cls)
        .join( ch_input_for_gsea.genesets )

    GSEA_GSEA(
        ch_input_for_gsea,
        ch_input_for_gsea.map{ index -> tuple(index[0].reference, index[0].target) },
        CUSTOM_TABULARTOGSEACHIP.out.chip.first()
    )
    // ----------------------------------------------------
    // Perform enrichment analysis with DECOUPLER
    // ----------------------------------------------------

    DECOUPLER_DECOUPLER(
        ch_input_for_other.input.filter{ index -> index[0].functional_method == 'decoupler' },
        ch_input_for_other.genesets.filter{ index -> index[0].functional_method == 'decoupler'},
        ch_input_for_other.features.filter{ index -> index[0].functional_method == 'decoupler'}
            .map{ meta, features_sheet, _features_id, _features_symbol -> [meta, features_sheet] }
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    PROPR_GREA(
        ch_input_for_other.input.filter{ index -> index[0].functional_method == 'grea' },
        ch_input_for_other.genesets.filter{ index -> index[0].functional_method == 'grea' }
    )

    // collect versions info
    ch_versions = ch_versions
        .mix(GPROFILER2_GOST.out.versions)
        .mix(CUSTOM_TABULARTOGSEAGCT.out.versions)
        .mix(CUSTOM_TABULARTOGSEACLS.out.versions)
        .mix(CUSTOM_TABULARTOGSEACHIP.out.versions)
        .mix(GSEA_GSEA.out.versions)
        .mix(DECOUPLER_DECOUPLER.out.versions)
        .mix(PROPR_GREA.out.versions)

    emit:
    // here we emit the outputs that will be useful afterwards in the
    // nf-core/differentialabundance pipeline

    // gprofiler2-specific outputs
    gprofiler2_all_enrich = GPROFILER2_GOST.out.all_enrich
    gprofiler2_sub_enrich = GPROFILER2_GOST.out.sub_enrich
    gprofiler2_plot_html  = GPROFILER2_GOST.out.plot_html

    // gsea-specific outputs
    gsea_report           = GSEA_GSEA.out.report_tsvs_ref.join(GSEA_GSEA.out.report_tsvs_target)

    // decoupler-specific outputs
    decoupler_dc_estimate = DECOUPLER_DECOUPLER.out.dc_estimate
    decoupler_dc_pvals = DECOUPLER_DECOUPLER.out.dc_pvals
    decoupler_png = DECOUPLER_DECOUPLER.out.png


    // grea-specific outputs
    grea_results          = PROPR_GREA.out.results

    // tool versions
    versions              = ch_versions
}
