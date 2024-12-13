
//
// Perform enrichment analysis
//
include { MYGENE } from "../../../modules/nf-core/mygene/main.nf"
include { PROPR_GREA as GREA } from "../../../modules/local/propr/grea/main.nf"
include { GPROFILER2_GOST } from "../../../modules/nf-core/gprofiler2/gost/main.nf"

include { GSEA_GSEA } from '../../../modules/nf-core/gsea/gsea/main.nf'
include { CUSTOM_TABULARTOGSEAGCT } from '../../../modules/nf-core/custom/tabulartogseagct/main.nf'
include { CUSTOM_TABULARTOGSEACLS } from '../../../modules/nf-core/custom/tabulartogseacls/main.nf'
include { TABULAR_TO_GSEA_CHIP } from '../../../modules/local/tabular_to_gsea_chip'

workflow ENRICHMENT {
    take:
    ch_counts                       // [ meta, counts] with meta keys: method, args_cor
    ch_results_genewise             // [ meta, results] with meta keys: method, args_cor
    ch_results_genewise_filtered    // [ meta, results] with meta keys: method, args_cor
    ch_adjacency                    // [ meta, adj_matrix] with meta keys: method, args_cor
    ch_contrasts
    ch_samplesheet
    ch_featuresheet
    ch_gene_sets

    // TODO: add ch_gm when provided by user, etc.

    main:

    // initialize empty results channels
    ch_enriched = Channel.empty()
    ch_gmt      = Channel.empty()
    ch_versions = Channel.empty()

    ch_counts
        .branch {
            gsea: it[0]["method"] == "gsea"
        }
        .set { ch_counts }

    ch_adjacency
        .branch {
            grea: it[0]["method"] == "grea"
        }
        .set { ch_adjacency }

    // ----------------------------------------------------
    // Construct gene set database
    // ----------------------------------------------------

    // TODO this should be optional, only run when there is no gene set data provided by user
    // TODO Use ch_gene_sets if it is not null
    MYGENE(ch_counts.take(1))  // only one data is provided to this pipeline
    ch_gmt      = MYGENE.out.gmt
    ch_versions = ch_versions.mix(MYGENE.out.versions)

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    GREA(ch_adjacency.grea.unique(), ch_gmt.collect())
    ch_enriched = ch_enriched.mix(GREA.out.results)
    ch_versions = ch_versions.mix(GREA.out.versions)

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // For GSEA, we need to convert normalised counts to a GCT format for
    // input, and process the sample sheet to generate class definitions
    // (CLS) for the variable used in each contrast

    CUSTOM_TABULARTOGSEAGCT ( ch_counts )

    // TODO: update CUSTOM_TABULARTOGSEACLS for value channel input per new
    // guidlines (rather than meta usage employed here)
    ch_contrasts_and_samples = ch_contrasts
        .map{it[0]} // revert back to contrasts meta map
        .combine( ch_samplesheet.map { it[1] } )

    CUSTOM_TABULARTOGSEACLS(ch_contrasts_and_samples)

    TABULAR_TO_GSEA_CHIP(
        ch_featuresheet.map{ it[1] },
        [params.features_id_col, params.features_name_col]
    )

    // The normalised matrix does not always have a contrast meta, so we
    // need a combine rather than a join here
    // Also add file name to metamap for easy access from modules.config
    // TODO combine the input channel with the ch_tools

    ch_gsea_inputs = CUSTOM_TABULARTOGSEAGCT.out.gct
        .map{ it.tail() }
        .combine(CUSTOM_TABULARTOGSEACLS.out.cls)
        .map{ tuple(it[1], it[0], it[2]) }
        .combine( ch_gmt.map { meta, gmt -> gmt } )

    println("__"+TABULAR_TO_GSEA_CHIP.out.chip)
    GSEA_GSEA(
        ch_gsea_inputs,
        ch_gsea_inputs.map{ tuple(it[0].reference, it[0].target) }, // *
        TABULAR_TO_GSEA_CHIP.out.chip.first()
    )

    // * Note: GSEA module currently uses a value channel for the mandatory
    // non-file arguments used to define contrasts, hence the indicated
    // usage of map to perform that transformation. An active subject of
    // debate
    GSEA_GSEA.out.report_tsvs_ref.view()
    ch_gsea_results = GSEA_GSEA.out.report_tsvs_ref
        .join(GSEA_GSEA.out.report_tsvs_target)

    ch_enriched = ch_enriched.combine(ch_gsea_results)


    // Record GSEA versions
    ch_versions = ch_versions
        .mix(TABULAR_TO_GSEA_CHIP.out.versions)
        .mix(GSEA_GSEA.out.versions)


    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    // todo: add gprofiler2 here

    // Define background file
    if (!params.gprofiler2_background_file) {
        // If deactivated, use empty list as "background"
        ch_background = []
    } else if (params.gprofiler2_background_file == "auto") {
        // If auto, use input matrix as background
        ch_background = ch_counts.map { meta, counts -> counts }
    } else {
        ch_background = Channel.from(file(params.gprofiler2_background_file, checkIfExists: true))
    }

    // rearrange channel for GPROFILER2_GOST process
    ch_gmt = ch_gmt.map { meta, gmt -> gmt }

    ch_results_genewise_filtered
        .branch {
            grea: it[0]["method"] == "gprofiler2"
        }
        .set { ch_results_genewise_filtered }

    GPROFILER2_GOST(ch_results_genewise_filtered, ch_gmt, ch_background)
    ch_versions = ch_versions.mix(GPROFILER2_GOST.out.versions)

    emit:
    enriched = ch_enriched  // channel: [ tsv ]
    versions = ch_versions  // channel: [ versions.yml ]
}
