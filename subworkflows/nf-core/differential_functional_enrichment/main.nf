
//
// Perform enrichment analysis
//
include { GPROFILER2_GOST          } from "../../../modules/nf-core/gprofiler2/gost/main.nf"
include { CUSTOM_TABULARTOGSEAGCT  } from '../../../modules/nf-core/custom/tabulartogseagct/main.nf'
include { CUSTOM_TABULARTOGSEACLS  } from '../../../modules/nf-core/custom/tabulartogseacls/main.nf'
include { CUSTOM_TABULARTOGSEACHIP } from '../../../modules/nf-core/custom/tabulartogseachip/main.nf'
include { GSEA_GSEA                } from '../../../modules/nf-core/gsea/gsea/main.nf'
include { PROPR_GREA               } from "../../../modules/nf-core/propr/grea/main.nf"

// Combine meta maps, including merging non-identical values of shared keys (e.g. 'id')
def mergeMaps(meta, meta2){
    (meta + meta2).collectEntries { k, v ->
        meta[k] && meta[k] != v ? [k, "${meta[k]}_${v}"] : [k, v]
    }
}

workflow DIFFERENTIAL_FUNCTIONAL_ENRICHMENT {
    take:
    // input data for functional analysis
    // They can be the results from differential expression analysis or abundance matrix
    // The functional analysis method to run should be explicitly provided
    ch_input                            // [ meta_input, input file, method to run ]

    // gene sets and background
    ch_gene_sets                        // [ meta_gmt, gmt file ]
    ch_background                       // [ meta_background, background file ]

    // other - for the moment these files are only needed for GSEA
    ch_contrasts                        // [ meta_contrast, contrast_variable, reference, target ]
    ch_samplesheet                      // [ meta_exp, samples sheet ]
    ch_featuresheet                     // [ meta_exp, features sheet, features id, features symbol ]

    main:

    ch_versions = Channel.empty()

    // Add method information into meta map of ch_input
    // This information is used later to determine which method to run for each input

    ch_input = ch_input.map {
        meta, file, analysis_method ->
        def meta_new = meta + [ 'method': analysis_method ]
        [ meta_new, file ]
    }

    // In the case of GSEA, it needs additional files coming from other channels that other methods don't use
    // here we define the input channel for the GSEA section

    def criteria = multiMapCriteria { meta_input, input, meta_exp, samplesheet, featuresheet, features_id, features_symbol, meta_contrasts, variable, reference, target ->
        def analysis_method = meta_input.method
        def meta_contrasts_new = meta_contrasts + [ 'variable': variable, 'reference': reference, 'target': target ]  // make sure variable, reference, target are in the meta
        def meta_all = mergeMaps(meta_contrasts_new, meta_input) + [ 'method': analysis_method ]
        input:
            [ meta_all, input ]
        contrasts_and_samples:
            [ meta_all, samplesheet ]
        features:
            [ meta_exp, featuresheet ]
        features_cols:
            [ features_id, features_symbol ]
    }
    ch_preinput_for_gsea = ch_input
        .filter { it[0].method == 'gsea' }
        .combine(ch_samplesheet.join(ch_featuresheet))
        .combine(ch_contrasts)
        .multiMap(criteria)

    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    GPROFILER2_GOST(
        ch_input.filter{ it[0].method == 'gprofiler2' },
        ch_gene_sets.collect(),
        ch_background.collect()
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // NOTE that GCT input can be more than 1, if they come from different tools (eg. limma, deseq2).
    // CLS input can be as many as combinations of input x contrasts
    // Whereas features can be only one file.

    CUSTOM_TABULARTOGSEAGCT(ch_preinput_for_gsea.input)

    CUSTOM_TABULARTOGSEACLS(ch_preinput_for_gsea.contrasts_and_samples)

    CUSTOM_TABULARTOGSEACHIP(
        ch_preinput_for_gsea.features.first(),
        ch_preinput_for_gsea.features_cols.first()
    )

    ch_input_for_gsea = CUSTOM_TABULARTOGSEAGCT.out.gct
        .join(CUSTOM_TABULARTOGSEACLS.out.cls)
        .combine( ch_gene_sets.map{it[1]}.collect() )

    GSEA_GSEA(
        ch_input_for_gsea,
        ch_input_for_gsea.map{ tuple(it[0].reference, it[0].target) },
        CUSTOM_TABULARTOGSEACHIP.out.chip.first()
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    PROPR_GREA(
        ch_input.filter{ it[0].method == 'grea' },
        ch_gene_sets.collect()
    )

    ch_input.view{ it -> "ch_input is ${it}"}
    ch_preinput_for_gsea.input.view { it -> "ch_preinputfor_gsea.input is ${it}"}
    ch_preinput_for_gsea.contrasts_and_samples.view { it -> "ch_preinput_for_gsea.contrasts_and_samples is ${it}"}
    ch_preinput_for_gsea.features.view{ it -> "ch_preinput_for_gsea.features is ${it}"}

    GPROFILER2_GOST.out.all_enrich.view { it -> "gprofiler2 output all enrich is ${it}"}
    GPROFILER2_GOST.out.sub_enrich.view { it -> "gprofiler2 output sub enrich is ${it}"}
    GPROFILER2_GOST.out.plot_html.view { it -> "gprofiler2 output html is ${it}"}

    CUSTOM_TABULARTOGSEAGCT.out.gct.view{ it -> "gct output is ${it}"}
    CUSTOM_TABULARTOGSEACLS.out.cls.view{ it -> "cls output is ${it}"}
    CUSTOM_TABULARTOGSEACHIP.out.chip.view{ it -> "chip output is ${it}"}
    GSEA_GSEA.out.report_tsvs_ref.view{ it -> "gsea ref report is ${it}"}
    GSEA_GSEA.out.report_tsvs_target.view{ it -> "gsea target report is ${it}"}

    PROPR_GREA.out.results.view{ it -> "grea output is ${it}"}

    ch_versions
        .mix(GPROFILER2_GOST.out.versions)
        .mix(CUSTOM_TABULARTOGSEAGCT.out.versions)
        .mix(CUSTOM_TABULARTOGSEACLS.out.versions)
        .mix(CUSTOM_TABULARTOGSEACHIP.out.versions)
        .mix(GSEA_GSEA.out.versions)
        .mix(PROPR_GREA.out.versions)
        .view{ it -> "versions is ${it}"}

    emit:
    // here we emit the outputs that will be useful afterwards in the
    // nf-core/differentialabundance pipeline

    // gprofiler2-specific outputs
    gprofiler2_all_enrich = GPROFILER2_GOST.out.all_enrich
    gprofiler2_sub_enrich = GPROFILER2_GOST.out.sub_enrich
    gprofiler2_plot_html  = GPROFILER2_GOST.out.plot_html

    // gsea-specific outputs
    gsea_report           = GSEA_GSEA.out.report_tsvs_ref
                                .join(GSEA_GSEA.out.report_tsvs_target)

    // grea-specific outputs
    grea_results          = PROPR_GREA.out.results

    // tool versions
    versions              = ch_versions
                                .mix(GPROFILER2_GOST.out.versions)
                                .mix(CUSTOM_TABULARTOGSEAGCT.out.versions)
                                .mix(CUSTOM_TABULARTOGSEACLS.out.versions)
                                .mix(CUSTOM_TABULARTOGSEACHIP.out.versions)
                                .mix(GSEA_GSEA.out.versions)
                                .mix(PROPR_GREA.out.versions)
}
