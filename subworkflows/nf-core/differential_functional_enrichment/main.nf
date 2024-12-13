
//
// Perform enrichment analysis
//
include { GPROFILER2_GOST         } from "../../../modules/nf-core/gprofiler2/gost/main.nf"
include { PROPR_GREA              } from "../../../modules/nf-core/propr/grea/main.nf"
// include { GSEA_GSEA               } from '../../../modules/nf-core/gsea/gsea/main.nf'
// include { CUSTOM_TABULARTOGSEAGCT } from '../../../modules/nf-core/custom/tabulartogseagct/main.nf'
// include { CUSTOM_TABULARTOGSEACLS } from '../../../modules/nf-core/custom/tabulartogseacls/main.nf'
// include { TABULAR_TO_GSEA_CHIP    } from '../../../modules/local/tabular_to_gsea_chip'

workflow ENRICHMENT {
    take:
    // input data for functional analysis
    // They can be the results from differential expression analysis or abundance matrix
    // The functional analysis method to run should be explicitly provided
    ch_input                            // [meta_input, input file, method to run]

    // gene sets and background
    ch_gene_sets                        // [ meta, gmt file]
    ch_background                       // [ meta, background file]

    main:

    // add the method information into the meta map
    // This information is used later to determine which method to run for each input

    ch_input = ch_input.map {
        meta, file, method ->
        def meta_new = meta + [ 'method': method ]
        [ meta_new, file ]
    }

    ch_input.view()

    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    GPROFILER2_GOST(
        ch_input.filter{ it[0].method == 'gprofiler2' }, 
        ch_gene_sets.map { meta, gmt -> gmt }, 
        ch_background
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    GREA(
        ch_input.filter{ it[0].method == 'grea' },
        ch_gene_sets.collect()
    )

    // // ----------------------------------------------------
    // // Perform enrichment analysis with GSEA
    // // ----------------------------------------------------

    // // For GSEA, we need to convert normalised counts to a GCT format for
    // // input, and process the sample sheet to generate class definitions
    // // (CLS) for the variable used in each contrast

    // CUSTOM_TABULARTOGSEAGCT ( ch_counts )

    // // TODO: update CUSTOM_TABULARTOGSEACLS for value channel input per new
    // // guidlines (rather than meta usage employed here)
    // ch_contrasts_and_samples = ch_contrasts
    //     .map{it[0]} // revert back to contrasts meta map
    //     .combine( ch_samplesheet.map { it[1] } )

    // CUSTOM_TABULARTOGSEACLS(ch_contrasts_and_samples)

    // TABULAR_TO_GSEA_CHIP(
    //     ch_featuresheet.map{ it[1] },
    //     [params.features_id_col, params.features_name_col]
    // )

    // // The normalised matrix does not always have a contrast meta, so we
    // // need a combine rather than a join here
    // // Also add file name to metamap for easy access from modules.config
    // // TODO combine the input channel with the ch_tools

    // ch_gsea_inputs = CUSTOM_TABULARTOGSEAGCT.out.gct
    //     .map{ it.tail() }
    //     .combine(CUSTOM_TABULARTOGSEACLS.out.cls)
    //     .map{ tuple(it[1], it[0], it[2]) }
    //     .combine( ch_gmt.map { meta, gmt -> gmt } )

    // println("__"+TABULAR_TO_GSEA_CHIP.out.chip)
    // GSEA_GSEA(
    //     ch_gsea_inputs,
    //     ch_gsea_inputs.map{ tuple(it[0].reference, it[0].target) }, // *
    //     TABULAR_TO_GSEA_CHIP.out.chip.first()
    // )

    // // * Note: GSEA module currently uses a value channel for the mandatory
    // // non-file arguments used to define contrasts, hence the indicated
    // // usage of map to perform that transformation. An active subject of
    // // debate
    // GSEA_GSEA.out.report_tsvs_ref.view()
    // ch_gsea_results = GSEA_GSEA.out.report_tsvs_ref
    //     .join(GSEA_GSEA.out.report_tsvs_target)

    // ch_enriched = ch_enriched.combine(ch_gsea_results)


    // // Record GSEA versions
    // ch_versions = ch_versions
    //     .mix(TABULAR_TO_GSEA_CHIP.out.versions)
    //     .mix(GSEA_GSEA.out.versions)

    // ----------------------------------------------------
    // Recollect results
    // ----------------------------------------------------

    // recollect enrichment main results from all tools

    ch_all_enrich =
        GPROFILER2_GOST.out.all_enrich
        .mix(GREA.out.results)

    ch_sub_enrich =
        GPROFILER2_GOST.out.sub_enrich

    // recollect plots needed for report

    ch_plot_html =
        GPROFILER2_GOST.out.plot_html

    // recollect versions

    ch_versions = GPROFILER2_GOST.out.versions
        .mix(GREA.out.versions)

    ch_all_enrich.view()

    emit:
    all_enrich = ch_all_enrich
    sub_enrich = ch_sub_enrich
    plot_html  = ch_plot_html
    versions   = ch_versions
}
