
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

    // add the method information into the meta map
    // This information is used later to determine which method to run for each input

    ch_input = ch_input.map {
        meta, file, analysis_method ->
        def meta_new = meta + [ 'method': analysis_method ]
        [ meta_new, file ]
    }

    // in the case of GSEA, it needs additional files coming from other channels that other methods don't use
    // here we define the input channel for the GSEA section

    def criteria = multiMapCriteria { meta_input, input, analysis_method, meta_exp, samplesheet, featuresheet, features_id, features_symbol, meta_contrasts, variable, reference, target ->
        input:
            meta_map = meta_input + [ 'method': analysis_method ]
            [ meta_map, input ]
        contrasts:
            meta_map = mergeMaps(meta_contrasts, meta_input) + [ 'method': analysis_method ]
            [ meta_map, variable, reference, target ]
        contrasts_and_samples:
            meta_map = mergeMaps(meta_contrasts, meta_exp) + [ 'method': analysis_method ]
            [ meta_map, samplesheet ]
        features:
            meta_map = meta_exp + [ 'method': analysis_method ]
            [ meta_map, featuresheet ]
        features_cols:
            [ features_id, features_symbol ]
    }

    ch_preinput_for_gsea = ch_input
        .filter { it[2] == 'gsea' }
        .combine(ch_samplesheet.join(ch_featuresheet))
        .combine(ch_contrasts)
        .multiMap(criteria)

    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    // deal with situations when some channels are empty
    if (ch_gene_sets == null) {
        ch_gene_sets_to_gprofiler2 = [[]]
    } else {
        ch_gene_sets_to_gprofiler2 = ch_gene_sets.map{ meta, gmt -> gmt }.collect()
    }
    if (ch_background == null) {
        ch_background_to_gprofiler2 = []
    } else {
        ch_background_to_gprofiler2 = ch_background.map{ meta, background -> background }.collect()
    }

    // TODO modify the module to take input files with meta values
    GPROFILER2_GOST(
        ch_input.filter{ it[0].method == 'gprofiler2' },
        ch_gene_sets_to_gprofiler2,
        ch_background_to_gprofiler2
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // For GSEA, we need to:
    //  - Convert normalised counts to a GCT format for input
    //  - Process the sample sheet to generate class definitions (CLS) for the variable used in each contrast
    //  - Process features sheet to generate CHIP file format

    // TODO: update CUSTOM_TABULARTOGSEACLS for value channel input per new
    // guidlines (rather than meta usage employed here)

    CUSTOM_TABULARTOGSEAGCT(ch_preinput_for_gsea.input)

    CUSTOM_TABULARTOGSEACLS(ch_preinput_for_gsea.contrasts_and_samples)

    CUSTOM_TABULARTOGSEACHIP(
        ch_preinput_for_gsea.features,
        ch_preinput_for_gsea.features_cols
    )

    if (ch_gene_sets == null) {
        ch_gene_sets_to_gsea = [[]]
    } else {
        ch_gene_sets_to_gsea = ch_gene_sets.map{ meta, gmt -> gmt }
    }

    ch_input_for_gsea = CUSTOM_TABULARTOGSEAGCT.out.gct
        .map{ it.tail() }
        .combine(CUSTOM_TABULARTOGSEACLS.out.cls)
        .map{ tuple(it[1], it[0], it[2]) }
        .combine( ch_gene_sets_to_gsea )

    GSEA_GSEA(
        ch_input_for_gsea,
        ch_input_for_gsea.map{ tuple(it[0].reference, it[0].target) },
        CUSTOM_TABULARTOGSEACHIP.out.chip.first()
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    if (ch_gene_sets == null) {
        ch_gene_sets_to_grea = [[], []]
    } else {
        ch_gene_sets_to_grea = ch_gene_sets.collect()
    }

    PROPR_GREA(
        ch_input.filter{ it[0].method == 'grea' },
        ch_gene_sets_to_grea
    )

    emit:
    // tool specific reports
    report_gprofiler2 = GPROFILER2_GOST.out.plot_html.map{it[1]}.flatMap().toList()
                            .combine(GPROFILER2_GOST.out.all_enrich.map{it[1]}.flatMap().toList())
                            .combine(GPROFILER2_GOST.out.sub_enrich.map{it[1]}.flatMap().toList())
    report_gsea       = GSEA_GSEA.out.report_tsvs_ref
                            .join(GSEA_GSEA.out.report_tsvs_target)
    report_grea       = PROPR_GREA.out.results

    // tool versions
    versions          = GPROFILER2_GOST.out.versions
                            .mix(CUSTOM_TABULARTOGSEAGCT.out.versions)
                            .mix(CUSTOM_TABULARTOGSEACLS.out.versions)
                            .mix(CUSTOM_TABULARTOGSEACHIP.out.versions)
                            .mix(GSEA_GSEA.out.versions)
                            .mix(PROPR_GREA.out.versions)
}
