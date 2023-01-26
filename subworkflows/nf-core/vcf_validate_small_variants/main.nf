include { RTGTOOLS_VCFEVAL      } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { RTGTOOLS_ROCPLOT      } from '../../../modules/nf-core/rtgtools/rocplot/main'
include { HAPPY_HAPPY           } from '../../../modules/nf-core/happy/happy/main'

workflow VCF_VALIDATE_SMALL_VARIANTS {

    take:
    ch_vcf                          // [mandatory] channel: [ meta, vcf, tbi, truth_vcf, truth_tbi, regions_bed, targets_bed ]
    ch_beds                         // [mandatory] channel: [ meta, regions_bed, targets_bed ]
    ch_fasta                        // [happy only] channel: [ meta, fasta ]
    ch_fasta_fai                    // [happy only] channel: [ meta, fasta_fai ]
    ch_vcfeval_sdf                  // [vcfeval only] channel: [ meta, sdf ]
    ch_happy_false_positive_regions // [optional] channel: [ meta, false_positives_bed ]
    ch_happy_stratification_tsv     // [optional] channel: [ meta, stratification_tsv ]
    ch_happy_stratification_beds    // [optional] channel: [ meta, [stratification_beds] ]
    tools                           // [mandatory] value: A comma-delimited list of the tools to use for validation (happy,vcfeval)

    main:

    ch_versions                             = Channel.empty()

    happy_vcf                               = Channel.empty()
    happy_tbi                               = Channel.empty()
    happy_indel_roc                         = Channel.empty()
    happy_indel_roc_pass                    = Channel.empty()
    happy_snp_roc                           = Channel.empty()
    happy_snp_roc_pass                      = Channel.empty()
    happy_roc                               = Channel.empty()
    happy_summary                           = Channel.empty()
    happy_extended_csv                      = Channel.empty()

    vcfeval_true_positive_vcf               = Channel.empty()
    vcfeval_true_positive_vcf_tbi           = Channel.empty()
    vcfeval_false_negative_vcf              = Channel.empty()
    vcfeval_false_negative_vcf_tbi          = Channel.empty()
    vcfeval_false_positive_vcf              = Channel.empty()
    vcfeval_false_positive_vcf_tbi          = Channel.empty()
    vcfeval_true_positive_baseline_vcf      = Channel.empty()
    vcfeval_true_positive_baseline_vcf_tbi  = Channel.empty()
    vcfeval_summary                         = Channel.empty()
    vcfeval_phasing                         = Channel.empty()
    vcfeval_snp_roc                         = Channel.empty()
    vcfeval_non_snp_roc                     = Channel.empty()
    vcfeval_weighted_roc                    = Channel.empty()

    rtgtools_snp_png_rocplot                = Channel.empty()
    rtgtools_non_snp_png_rocplot            = Channel.empty()
    rtgtools_weighted_png_rocplot           = Channel.empty()

    rtgtools_snp_svg_rocplot                = Channel.empty()
    rtgtools_non_snp_svg_rocplot            = Channel.empty()
    rtgtools_weighted_svg_rocplot           = Channel.empty()

    list_tools = tools.tokenize(",")

    ch_input = ch_vcf.join(ch_beds)

    if("happy" in list_tools){
        happy_input = ch_input
            .map { meta, vcf, tbi, truth_vcf, truth_tbi, regions_bed, targets_bed ->
                [ meta, vcf, truth_vcf, regions_bed, targets_bed ]
            }
        HAPPY_HAPPY (
            happy_input,
            ch_fasta,
            ch_fasta_fai,
            ch_happy_false_positive_regions,
            ch_happy_stratification_tsv,
            ch_happy_stratification_beds
        )
        ch_versions = ch_versions.mix(HAPPY_HAPPY.out.versions.first())

        happy_vcf               = HAPPY_HAPPY.out.vcf
        happy_tbi               = HAPPY_HAPPY.out.tbi
        happy_indel_roc         = HAPPY_HAPPY.out.roc_indel_locations_csv
        happy_indel_roc_pass    = HAPPY_HAPPY.out.roc_indel_locations_pass_csv
        happy_snp_roc           = HAPPY_HAPPY.out.roc_snp_locations_csv
        happy_snp_roc_pass      = HAPPY_HAPPY.out.roc_snp_locations_pass_csv
        happy_roc               = HAPPY_HAPPY.out.roc_all_csv
        happy_summary           = HAPPY_HAPPY.out.summary_csv
        happy_extended_csv      = HAPPY_HAPPY.out.extended_csv
    }

    if("vcfeval" in list_tools){
        RTGTOOLS_VCFEVAL(
            ch_input.map { it[0..-2] + [[]] },
            ch_vcfeval_sdf
        )
        ch_versions = ch_versions.mix(RTGTOOLS_VCFEVAL.out.versions.first())

        ch_rocplot_input = RTGTOOLS_VCFEVAL.out.snp_roc
            .map { meta, tsv ->
                [ meta + [roc_type:'snp'], tsv ]
            }
            .mix(
                RTGTOOLS_VCFEVAL.out.non_snp_roc.map { meta, tsv ->
                    [ meta + [roc_type:'non_snp'], tsv ]
                },
                RTGTOOLS_VCFEVAL.out.weighted_roc.map { meta, tsv ->
                    [ meta + [roc_type:'weighted'], tsv ]
                }
            )

        vcfeval_true_positive_vcf               = RTGTOOLS_VCFEVAL.out.tp_vcf
        vcfeval_true_positive_vcf_tbi           = RTGTOOLS_VCFEVAL.out.tp_tbi
        vcfeval_false_negative_vcf              = RTGTOOLS_VCFEVAL.out.fn_vcf
        vcfeval_false_negative_vcf_tbi          = RTGTOOLS_VCFEVAL.out.fn_tbi
        vcfeval_false_positive_vcf              = RTGTOOLS_VCFEVAL.out.fp_vcf
        vcfeval_false_positive_vcf_tbi          = RTGTOOLS_VCFEVAL.out.fp_tbi
        vcfeval_true_positive_baseline_vcf      = RTGTOOLS_VCFEVAL.out.baseline_vcf
        vcfeval_true_positive_baseline_vcf_tbi  = RTGTOOLS_VCFEVAL.out.baseline_tbi
        vcfeval_summary                         = RTGTOOLS_VCFEVAL.out.summary
        vcfeval_phasing                         = RTGTOOLS_VCFEVAL.out.phasing
        vcfeval_snp_roc                         = RTGTOOLS_VCFEVAL.out.snp_roc
        vcfeval_non_snp_roc                     = RTGTOOLS_VCFEVAL.out.non_snp_roc
        vcfeval_weighted_roc                    = RTGTOOLS_VCFEVAL.out.weighted_roc

        RTGTOOLS_ROCPLOT(
            ch_rocplot_input
        )

        ch_versions = ch_versions.mix(RTGTOOLS_ROCPLOT.out.versions.first())

        rocplot_out_png = RTGTOOLS_ROCPLOT.out.png
            .branch { meta, png ->
                roc_type = meta.roc_type
                new_meta = meta.findAll { !(it.key == "roc_type") }

                snp:        roc_type == "snp"
                non_snp:    roc_type == "non_snp"
                weighted:   roc_type == "weighted"
            }

        rocplot_out_svg = RTGTOOLS_ROCPLOT.out.svg
            .branch { meta, svg ->
                roc_type = meta.roc_type
                new_meta = meta.findAll { !(it.key == "roc_type") }

                snp:        roc_type == "snp"
                non_snp:    roc_type == "non_snp"
                weighted:   roc_type == "weighted"
            }

        rtgtools_snp_png_rocplot        = rocplot_out_png.snp
        rtgtools_non_snp_png_rocplot    = rocplot_out_png.non_snp
        rtgtools_weighted_png_rocplot   = rocplot_out_png.weighted

        rtgtools_snp_svg_rocplot        = rocplot_out_svg.snp
        rtgtools_non_snp_svg_rocplot    = rocplot_out_svg.non_snp
        rtgtools_weighted_svg_rocplot   = rocplot_out_svg.weighted

    }

    emit:
    happy_vcf                               // channel: [ meta, vcf ]
    happy_tbi                               // channel: [ meta, tbi ]
    happy_indel_roc                         // channel: [ meta, csv ]
    happy_indel_roc_pass                    // channel: [ meta, csv ]
    happy_snp_roc                           // channel: [ meta, csv ]
    happy_snp_roc_pass                      // channel: [ meta, csv ]
    happy_roc                               // channel: [ meta, csv ]
    happy_summary                           // channel: [ meta, csv ]
    happy_extended_csv                      // channel: [ meta, csv ]

    vcfeval_true_positive_vcf               // channel: [ meta, vcf ]
    vcfeval_true_positive_vcf_tbi           // channel: [ meta, tbi ]
    vcfeval_false_negative_vcf              // channel: [ meta, vcf ]
    vcfeval_false_negative_vcf_tbi          // channel: [ meta, tbi ]
    vcfeval_false_positive_vcf              // channel: [ meta, vcf ]
    vcfeval_false_positive_vcf_tbi          // channel: [ meta, tbi ]
    vcfeval_true_positive_baseline_vcf      // channel: [ meta, vcf ]
    vcfeval_true_positive_baseline_vcf_tbi  // channel: [ meta, tbi ]
    vcfeval_summary                         // channel: [ meta, summary ]
    vcfeval_phasing                         // channel: [ meta, phasing ]
    vcfeval_snp_roc                         // channel: [ meta, tsv ]
    vcfeval_non_snp_roc                     // channel: [ meta, tsv ]
    vcfeval_weighted_roc                    // channel: [ meta, tsv ]

    rtgtools_snp_png_rocplot                // channel: [ meta, png ]
    rtgtools_non_snp_png_rocplot            // channel: [ meta, png ]
    rtgtools_weighted_png_rocplot           // channel: [ meta, png ]
    rtgtools_snp_svg_rocplot                // channel: [ meta, svg ]
    rtgtools_non_snp_svg_rocplot            // channel: [ meta, svg ]
    rtgtools_weighted_svg_rocplot           // channel: [ meta, svg ]

    versions = ch_versions                  // channel: [ versions.yml ]
}

