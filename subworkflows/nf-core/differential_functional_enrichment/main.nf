
//
// Perform enrichment analysis
//
include { GPROFILER2_GOST          } from "../../../modules/nf-core/gprofiler2/gost/main.nf"
include { CUSTOM_TABULARTOGSEAGCT  } from '../../../modules/nf-core/custom/tabulartogseagct/main.nf'
include { CUSTOM_TABULARTOGSEACLS  } from '../../../modules/nf-core/custom/tabulartogseacls/main.nf'
include { CUSTOM_TABULARTOGSEACHIP } from '../../../modules/nf-core/custom/tabulartogseachip/main.nf'
include { GSEA_GSEA                } from '../../../modules/nf-core/gsea/gsea/main.nf'
include { PROPR_GREA               } from "../../../modules/nf-core/propr/grea/main.nf"

workflow DIFFERENTIAL_FUNCTIONAL_ENRICHMENT {
    take:
    // input data for functional analysis
    // They can be the results from differential expression analysis or abundance matrix
    // The functional analysis method to run should be explicitly provided
    ch_input                            // [ meta_input, input file, method to run ]

    // gene sets and background
    ch_gene_sets                        // [ meta, gmt file ]
    ch_background                       // [ background file ]

    // other
    ch_contrasts                        // [ meta_contrast, contrast_variable, reference, target ]
    ch_samplesheet                      // [ meta, samples sheet ]
    ch_featuresheet                     // [ meta, features sheet, features id, features symbol ]

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

    // TODO modify these modules to take input files with meta values

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

    GPROFILER2_GOST(
        ch_input.filter{ it[0].method == 'gprofiler2' },
        ch_gene_sets_to_gprofiler2,
        ch_background_to_gprofiler2
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    ch_report_gsea = Channel.empty()
    ch_versions_gsea = Channel.empty()

    // count elements that need to run through gsea
    // if zero, skip
    def count_gsea = 0
    ch_input.filter{ it[0].method == 'gsea' }
        .view { item ->
            count_gsea++
            return null  // Don't print anything
        }

    if (count_gsea > 0) {

        // For GSEA, we need to:
        //  - Convert normalised counts to a GCT format for input
        //  - Process the sample sheet to generate class definitions (CLS) for the variable used in each contrast
        //  - Process features sheet to generate CHIP file format

        // TODO: update CUSTOM_TABULARTOGSEACLS for value channel input per new
        // guidlines (rather than meta usage employed here)

        CUSTOM_TABULARTOGSEAGCT( ch_input.filter{ it[0].method == 'grea' } )

        ch_contrasts_and_samples = ch_contrasts
            .map{ meta_contrast, contrast_variable, reference, target ->
                meta_contrast['variable'] = contrast_variable  // make sure meta contains the contrast variable
                return meta_contrast
            }
            .combine( ch_samplesheet.map { meta, samplesheet -> samplesheet } )
        CUSTOM_TABULARTOGSEACLS(ch_contrasts_and_samples)

        ch_features_to_chip = ch_featuresheet
            .multiMap { meta, features, features_id_col, features_name_col ->
                tabular: [meta, features]
                cols: [features_id_col, features_name_col]
            }
        CUSTOM_TABULARTOGSEACHIP(
            ch_features_to_chip.tabular,
            ch_features_to_chip.cols
        )

        // The normalised matrix does not always have a contrast meta, so we
        // need a combine rather than a join here
        // Also add file name to metamap for easy access from modules.config

        ch_gsea_inputs = CUSTOM_TABULARTOGSEAGCT.out.gct
            .map{ it.tail() }
            .combine(CUSTOM_TABULARTOGSEACLS.out.cls)
            .map{ tuple(it[1], it[0], it[2]) }
            .combine( ch_gene_sets.map{ meta, gmt -> gmt } )

        println("__"+CUSTOM_TABULARTOGSEACHIP.out.chip)
        GSEA_GSEA(
            ch_gsea_inputs,
            ch_gsea_inputs.map{ tuple(it[0].reference, it[0].target) }, // *
            CUSTOM_TABULARTOGSEACHIP.out.chip.first()
        )

        ch_report_gsea = GSEA_GSEA.out.report_tsvs_ref
                            .join(GSEA_GSEA.out.report_tsvs_target)
        ch_versions_gsea = CUSTOM_TABULARTOGSEAGCT.out.versions
                            .mix(CUSTOM_TABULARTOGSEACLS.out.versions)
                            .mix(CUSTOM_TABULARTOGSEACHIP.out.versions)
                            .mix(GSEA_GSEA.out.versions)
    }

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
    // report_gsea       = ch_report_gsea
    report_grea       = PROPR_GREA.out.results

    // tool versions
    versions          = GPROFILER2_GOST.out.versions
                            // .mix(ch_versions_gsea)
                            .mix(PROPR_GREA.out.versions)
}
