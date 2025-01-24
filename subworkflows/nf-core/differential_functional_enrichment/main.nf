
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
    ch_input                 // [ meta_input, input file, genesets file, background file, method to run ]

    // other - for the moment these files are only needed for GSEA
    // as it is the only one that takes expression data as input
    // if in the future this setting is changed, this section could be removed
    ch_contrasts             // [ meta_contrast, contrast_variable, reference, target ]
    ch_samplesheet           // [ meta_exp, samples sheet ]
    ch_featuresheet          // [ meta_exp, features sheet, features id, features symbol ]

    main:

    ch_versions = Channel.empty()

    // Add method information into meta map of ch_input
    // This information is used later to determine which method to run for each input
    // Also, reorganize the structure to match them with the modules' input organization

    ch_input = ch_input
        .multiMap {
            meta_input, file, genesets, background, analysis_method ->
            def meta_new = meta_input + [ 'method': analysis_method ]
            input:
                [ meta_new, file ]
            genesets:
                [ meta_new, genesets ]  // NOTE here we assume that the modules will not make use of meta_genesets and meta_background
            background:
                [ meta_new, background ]
        }

    // In the case of GSEA, it needs additional files coming from other channels that other methods don't use
    // here we define the input channel for the GSEA section

    def criteria = multiMapCriteria { meta_input, input, genesets, meta_exp, samplesheet, featuresheet, features_id, features_symbol, meta_contrasts, variable, reference, target ->
        def meta_contrasts_new = meta_contrasts + [ 'variable': variable, 'reference': reference, 'target': target ]  // make sure variable, reference, target are in the meta
        def meta_all = mergeMaps(meta_contrasts_new, meta_input)
        input:
            [ meta_all, input ]
        genesets:
            [ meta_all, genesets ]
        contrasts_and_samples:
            [ meta_all, samplesheet ]
        features:
            [ meta_exp, featuresheet ]
        features_cols:
            [ features_id, features_symbol ]
    }
    ch_preinput_for_gsea = ch_input.input
        .join(ch_input.genesets)
        .filter{ it[0].method == 'gsea' }
        .combine(ch_samplesheet.join(ch_featuresheet))
        .combine(ch_contrasts)
        .multiMap(criteria)

    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    GPROFILER2_GOST(
        ch_input.input.filter{ it[0].method == 'gprofiler2' },
        ch_input.genesets.filter{ it[0].method == 'gprofiler2'},
        ch_input.background.filter{ it[0].method == 'gprofiler2'}
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
        .join( ch_preinput_for_gsea.genesets )

    GSEA_GSEA(
        ch_input_for_gsea,
        ch_input_for_gsea.map{ tuple(it[0].reference, it[0].target) },
        CUSTOM_TABULARTOGSEACHIP.out.chip.first()
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    PROPR_GREA(
        ch_input.input.filter{ it[0].method == 'grea' },
        ch_input.genesets.filter{ it[0].method == 'grea' }
    )

    // collect versions info
    ch_versions = ch_versions
        .mix(GPROFILER2_GOST.out.versions)
        .mix(CUSTOM_TABULARTOGSEAGCT.out.versions)
        .mix(CUSTOM_TABULARTOGSEACLS.out.versions)
        .mix(CUSTOM_TABULARTOGSEACHIP.out.versions)
        .mix(GSEA_GSEA.out.versions)
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

    // grea-specific outputs
    grea_results          = PROPR_GREA.out.results

    // tool versions
    versions              = ch_versions
}
