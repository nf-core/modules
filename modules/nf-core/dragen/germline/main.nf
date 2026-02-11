process DRAGEN {
    tag "$meta.id"
    label 'process_dragen'

    // ATTENTION: No conda env or container image as Dragen requires specialized hardware to run

    input:
    tuple val(meta), path(input)
    path checkfingerprint_expected_vcf
    path cnv_combined_counts
    path cnv_exclude_bed
    path cnv_population_b_allele_vcf
    path cnv_segmentation_bed
    path cnv_target_bed
    path cram_reference
    path dbsnp
    path fastqc_adapter_file
    path fastqc_kmer_file
    path ora_reference
    path qc_coverage_region
    path qc_cross_cont_vcf
    path ref_dir
    path repeat_genotype_ref_fasta
    path repeat_genotype_specs
    path sv_call_regions_bed
    path sv_exclusion_bed
    path sv_forcegt_vcf
    path sv_systematic_noise
    path trim_adapter_read
    path trim_adapter_read_5prime
    path variant_annotation_data
    path vc_combine_phased_variants_distance_bed
    path vc_excluded_regions_bed
    path vc_forcegt_vcf
    path vc_log_bed
    path vc_mapping_metrics
    path vc_ml_dir
    path vc_ntd_error_params
    path vc_roh_blacklist_bed
    path vc_snp_error_cal_bed
    path vc_systematic_noise
    path vc_target_bed
    path vd_eh_vcf
    path vd_small_variant_vcf
    path vd_sv_vcf
    path vntr_catalog_bed

    output:
    tuple val(meta), path("${prefix}-replay.json")                           , optional: true, emit: replay_json
    tuple val(meta), path("${prefix}.*_coverage_metrics.csv")                , optional: true, emit: coverage_metrics_csv
    tuple val(meta), path("${prefix}.SJ.saturation.txt")                     , optional: true, emit: SJ_saturation_txt
    tuple val(meta), path("${prefix}.baf.bedgraph.gz")                       , optional: true, emit: baf_bedgraph_gz
    tuple val(meta), path("${prefix}.baf.bw")                                , optional: true, emit: baf_bw
    tuple val(meta), path("${prefix}.baf.seg")                               , optional: true, emit: baf_seg
    tuple val(meta), path("${prefix}.baf.seq.bw")                            , optional: true, emit: baf_seq_bw
    tuple val(meta), path("${prefix}.ballele.counts.gz")                     , optional: true, emit: ballele_counts_gz
    tuple val(meta), path("${prefix}.bam")                                   , optional: true, emit: bam
    tuple val(meta), path("${prefix}.bam.bai")                               , optional: true, emit: bam_bai
    tuple val(meta), path("${prefix}.bam.md5sum")                            , optional: true, emit: bam_md5sum
    tuple val(meta), path("${prefix}.cnv.excluded_intervals.bed.gz")         , optional: true, emit: cnv_excluded_intervals_bed_gz
    tuple val(meta), path("${prefix}.cnv.gff3")                              , optional: true, emit: cnv_gff3
    tuple val(meta), path("${prefix}.cnv.igv_session.xml")                   , optional: true, emit: cnv_igv_session_xml
    tuple val(meta), path("${prefix}.cnv.pon_correlation.txt.gz")            , optional: true, emit: cnv_pon_correlation_txt_gz
    tuple val(meta), path("${prefix}.cnv.pon_metrics.tsv.gz")                , optional: true, emit: cnv_pon_metrics_tsv_gz
    tuple val(meta), path("${prefix}.cnv.purity.coverage.models.tsv")        , optional: true, emit: cnv_purity_coverage_models_tsv
    tuple val(meta), path("${prefix}.cnv.segdups.joint_coverage.tsv.gz")     , optional: true, emit: cnv_segdups_joint_coverage_tsv_gz
    tuple val(meta), path("${prefix}.cnv.segdups.rescued_intervals.tsv.gz")  , optional: true, emit: cnv_segdups_rescued_intervals_tsv_gz
    tuple val(meta), path("${prefix}.cnv.segdups.site_ratios.tsv.gz")        , optional: true, emit: cnv_segdups_site_ratios_tsv_gz
    tuple val(meta), path("${prefix}.cnv.vcf")                               , optional: true, emit: cnv_vcf
    tuple val(meta), path("${prefix}.cnv.vcf.gz")                            , optional: true, emit: cnv_vcf_gz
    tuple val(meta), path("${prefix}.cnv.vcf.gz.md5sum")                     , optional: true, emit: cnv_vcf_gz_md5sum
    tuple val(meta), path("${prefix}.cnv.vcf.gz.tbi")                        , optional: true, emit: cnv_vcf_gz_tbi
    tuple val(meta), path("${prefix}.cnv_metrics.csv")                       , optional: true, emit: cnv_metrics_csv
    tuple val(meta), path("${prefix}.cnv_sv.vcf.gz")                         , optional: true, emit: cnv_sv_vcf_gz
    tuple val(meta), path("${prefix}.combined.counts.txt.gz")                , optional: true, emit: combined_counts_txt_gz
    tuple val(meta), path("${prefix}.cram")                                  , optional: true, emit: cram
    tuple val(meta), path("${prefix}.cram.crai")                             , optional: true, emit: cram_crai
    tuple val(meta), path("${prefix}.cram.md5sum")                           , optional: true, emit: cram_md5sum
    tuple val(meta), path("${prefix}.excluded_intervals.bed.gz")             , optional: true, emit: excluded_intervals_bed_gz
    tuple val(meta), path("${prefix}.fastqc_metrics.csv")                    , optional: true, emit: fastqc_metrics_csv
    tuple val(meta), path("${prefix}.fragment_length_hist.csv")              , optional: true, emit: fragment_length_hist_csv
    tuple val(meta), path("${prefix}.g.vcf.gz")                              , optional: true, emit: g_vcf_gz
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi")                          , optional: true, emit: g_vcf_gz_tbi
    tuple val(meta), path("${prefix}.gc_metrics.csv")                        , optional: true, emit: gc_metrics_csv
    tuple val(meta), path("${prefix}.hard-filtered.baf.bw")                  , optional: true, emit: hard_filtered_baf_bw
    tuple val(meta), path("${prefix}.hard-filtered.gvcf.gz.md5sum")          , optional: true, emit: hard_filtered_gvcf_gz_md5sum
    tuple val(meta), path("${prefix}.hard-filtered.vcf")                     , optional: true, emit: hard_filtered_vcf
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz")                  , optional: true, emit: hard_filtered_vcf_gz
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz.md5sum")           , optional: true, emit: hard_filtered_vcf_gz_md5sum
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz.tbi")              , optional: true, emit: hard_filtered_vcf_gz_tbi
    tuple val(meta), path("${prefix}.improper.pairs.bw")                     , optional: true, emit: improper_pairs_bw
    tuple val(meta), path("${prefix}.impute.chunk.out.txt")                  , optional: true, emit: impute_chunk_out_txt
    tuple val(meta), path("${prefix}.impute.phase.out.txt")                  , optional: true, emit: impute_phase_out_txt
    tuple val(meta), path("${prefix}.impute.vcf.gz")                         , optional: true, emit: impute_vcf_gz
    tuple val(meta), path("${prefix}.insert-stats.tab")                      , optional: true, emit: insert_stats_tab
    tuple val(meta), path("${prefix}.mapping_metrics.csv")                   , optional: true, emit: mapping_metrics_csv
    tuple val(meta), path("${prefix}.multiomics.barcodeSummary.tsv")         , optional: true, emit: multiomics_barcodeSummary_tsv
    tuple val(meta), path("${prefix}.multiomics.barcodes.tsv.gz")            , optional: true, emit: multiomics_barcodes_tsv_gz
    tuple val(meta), path("${prefix}.multiomics.features.tsv.gz")            , optional: true, emit: multiomics_features_tsv_gz
    tuple val(meta), path("${prefix}.multiomics.filtered.barcodes.tsv.gz")   , optional: true, emit: multiomics_filtered_barcodes_tsv_gz
    tuple val(meta), path("${prefix}.multiomics.matrix.mtx.gz")              , optional: true, emit: multiomics_matrix_mtx_gz
    tuple val(meta), path("${prefix}.multiomics.metrics.csv")                , optional: true, emit: multiomics_metrics_csv
    tuple val(meta), path("${prefix}.pcr-mode*.log")                         , optional: true, emit: pcr_model_log
    tuple val(meta), path("${prefix}.ploidy.vcf")                            , optional: true, emit: ploidy_vcf
    tuple val(meta), path("${prefix}.ploidy.vcf.gz")                         , optional: true, emit: ploidy_vcf_gz
    tuple val(meta), path("${prefix}.ploidy.vcf.gz.md5sum")                  , optional: true, emit: ploidy_vcf_gz_md5sum
    tuple val(meta), path("${prefix}.ploidy.vcf.gz.tbi")                     , optional: true, emit: ploidy_vcf_gz_tbi
    tuple val(meta), path("${prefix}.ploidy_estimation_metrics.csv")         , optional: true, emit: ploidy_estimation_metrics_csv
    tuple val(meta), path("${prefix}.preprocess.vcf.gz")                     , optional: true, emit: preprocess_vcf_gz
    tuple val(meta), path("${prefix}.quant.genes.sf")                        , optional: true, emit: quant_genes_sf
    tuple val(meta), path("${prefix}.quant.metrics.csv")                     , optional: true, emit: quant_metrics_csv
    tuple val(meta), path("${prefix}.quant.sf")                              , optional: true, emit: quant_sf
    tuple val(meta), path("${prefix}.quant.transcript_coverage.txt")         , optional: true, emit: quant_transcript_coverage_txt
    tuple val(meta), path("${prefix}.quant.transcript_fragment_lengths.txt") , optional: true, emit: quant_transcript_fragment_lengths_txt
    tuple val(meta), path("${prefix}.realigned-regions.bed")                 , optional: true, emit: realigned_regions_bed
    tuple val(meta), path("${prefix}.repeats.bam")                           , optional: true, emit: repeats_bam
    tuple val(meta), path("${prefix}.repeats.vcf")                           , optional: true, emit: repeats_vcf
    tuple val(meta), path("${prefix}.repeats.vcf.gz")                        , optional: true, emit: repeats_vcf_gz
    tuple val(meta), path("${prefix}.repeats.vcf.gz.tbi")                    , optional: true, emit: repeats_vcf_gz_tbi
    tuple val(meta), path("${prefix}.roh.bed")                               , optional: true, emit: roh_bed
    tuple val(meta), path("${prefix}.roh_metrics.csv")                       , optional: true, emit: roh_metrics_csv
    tuple val(meta), path("${prefix}.sam")                                   , optional: true, emit: sam
    tuple val(meta), path("${prefix}.seg")                                   , optional: true, emit: seg
    tuple val(meta), path("${prefix}.seg.bw")                                , optional: true, emit: seg_bw
    tuple val(meta), path("${prefix}.seg.c.partial.Y")                       , optional: true, emit: seg_c_partial_y
    tuple val(meta), path("${prefix}.seg.called")                            , optional: true, emit: seg_called
    tuple val(meta), path("${prefix}.seg.called.merged")                     , optional: true, emit: seg_called_merged
    tuple val(meta), path("${prefix}.select.gvcf")                           , optional: true, emit: select_gvcf
    tuple val(meta), path("${prefix}.small_indel_dedup")                     , optional: true, emit: small_indel_dedup
    tuple val(meta), path("${prefix}.smn_dedup")                             , optional: true, emit: smn_dedup
    tuple val(meta), path("${prefix}.snperror-sampler.log")                  , optional: true, emit: snperror_samples_log
    tuple val(meta), path("${prefix}.sv.vcf")                                , optional: true, emit: sv_vcf
    tuple val(meta), path("${prefix}.sv.vcf.gz.tbi")                         , optional: true, emit: sv_vcf_gz_tbi
    tuple val(meta), path("${prefix}.sv_metrics.csv")                        , optional: true, emit: sv_metrics_csv
    tuple val(meta), path("${prefix}.target.counts.bw")                      , optional: true, emit: target_counts_bw
    tuple val(meta), path("${prefix}.target.counts.diploid.bw")              , optional: true, emit: target_counts_diploid_bw
    tuple val(meta), path("${prefix}.target.counts.gc-corrected.gz")         , optional: true, emit: target_counts_gc_corrected_gz
    tuple val(meta), path("${prefix}.target.counts.gz")                      , optional: true, emit: target_counts_gz
    tuple val(meta), path("${prefix}.targeted.json")                         , optional: true, emit: targeted_json
    tuple val(meta), path("${prefix}.targeted.vcf.gz")                       , optional: true, emit: targeted_vcf_gz
    tuple val(meta), path("${prefix}.time_metrics.csv")                      , optional: true, emit: time_metrics_csv
    tuple val(meta), path("${prefix}.tn.bw")                                 , optional: true, emit: tn_bw
    tuple val(meta), path("${prefix}.tn.tsv.gz")                             , optional: true, emit: tn_tsv_gz
    tuple val(meta), path("${prefix}.trimmer_metrics.csv")                   , optional: true, emit: trimmer_metrics_csv
    tuple val(meta), path("${prefix}.vc_hethom_ratio_metrics.csv")           , optional: true, emit: vc_hethom_ratio_metrics_csv
    tuple val(meta), path("${prefix}.vc_metrics.csv")                        , optional: true, emit: vc_metrics_csv
    tuple val(meta), path("${prefix}.vcf.gz")                                , optional: true, emit: vcf_gz
    tuple val(meta), path("${prefix}.vcf.gz.md5sum")                         , optional: true, emit: vcf_gz_md5sum
    tuple val(meta), path("${prefix}.vcf.gz.tbi")                            , optional: true, emit: vcf_gz_tbi
    tuple val(meta), path("${prefix}.wgs_contig_mean_cov.csv")               , optional: true, emit: wgs_contig_mean_cov_csv
    tuple val(meta), path("${prefix}.wgs_coverage_metrics.csv")              , optional: true, emit: wgs_coverage_metrics_csv
    tuple val(meta), path("${prefix}.wgs_fine_hist.csv")                     , optional: true, emit: wgs_fine_hist_csv
    tuple val(meta), path("${prefix}.wgs_hist.csv")                          , optional: true, emit: wgs_hist_csv
    tuple val(meta), path("${prefix}.wgs_overall_mean_cov.csv")              , optional: true, emit: wgs_overall_mean_cov_csv
    tuple val(meta), path("${prefix}_callability.bed")                       , optional: true, emit: callability_bed
    tuple val(meta), path("${prefix}_chr_start-end.impute.phase.vcf.gz")     , optional: true, emit: chr_start_end_impute_phase_vcf_gz
    tuple val(meta), path("${prefix}_contig_mean_cov.csv")                   , optional: true, emit: contig_mean_cov_csv
    tuple val(meta), path("${prefix}_cov_report.bed")                        , optional: true, emit: cov_report_bed
    tuple val(meta), path("${prefix}_evidence.{b,cr,s}am")                   , optional: true, emit: evidence
    tuple val(meta), path("${prefix}_fine_hist.csv")                         , optional: true, emit: fine_hist_csv
    tuple val(meta), path("${prefix}_full_res.bed")                          , optional: true, emit: full_res_bed
    tuple val(meta), path("${prefix}_hist.csv")                              , optional: true, emit: hist_csv
    tuple val(meta), path("${prefix}_overall_mean_cov.csv")                  , optional: true, emit: overall_mean_cov_csv
    tuple val(meta), path("${prefix}_read_cov_report.bed")                   , optional: true, emit: read_cov_report_bed
    path "./sort_spill/partitions.txt"                                       , optional: true, emit: partitions_txt
    path "./sv/results/stats/alignmentStatsSummary.txt"                      , optional: true, emit: alignmentStatsSummary_txt
    path "./sv/results/stats/candidate_metrics.csv"                          , optional: true, emit: candidate_metrics_csv
    path "./sv/results/stats/diploidSV.sv_metrics.csv"                       , optional: true, emit: diploidSV_sv_metrics_csv
    path "./sv/results/stats/graph_metrics.csv"                              , optional: true, emit: graph_metrics_csv
    path "./sv/results/stats/svCandidateGenerationStats.tsv"                 , optional: true, emit: svCandidateGenerationStats_tsv
    path "./sv/results/stats/svCandidateGenerationStats.xml"                 , optional: true, emit: svCandidateGenerationStats_xml
    path "./sv/results/stats/svLocusGraphStats.tsv"                          , optional: true, emit: svLocusGraphStats_tsv
    path "./sv/results/variants/candidateSV.vcf.gz"                          , optional: true, emit: candidateSV_vcf_gz
    path "./sv/results/variants/candidateSV.vcf.gz.tbi"                      , optional: true, emit: candidateSV_vcf_gz_tbi
    path "./sv/results/variants/diploidSV.vcf.gz"                            , optional: true, emit: diploidSV_vcf_gz
    path "./sv/results/variants/diploidSV.vcf.gz.tbi"                        , optional: true, emit: diploidSV_vcf_gz_tbi
    path "./sv/workspace/alignmentStats.xml"                                 , optional: true, emit: alignmentStats_xml
    path "./sv/workspace/chromDepth.txt"                                     , optional: true, emit: chromDepth_txt
    path "./sv/workspace/edgeRuntimeLog.txt"                                 , optional: true, emit: edgeRuntimeLog_txt
    path "./sv/workspace/genomeSegmentScanDebugInfo.txt"                     , optional: true, emit: genomeSegmentScanDebugInfo_txt
    path "./sv/workspace/logs/config_log.txt"                                , optional: true, emit: config_log_txt
    path "./sv/workspace/svLocusGraph.bin"                                   , optional: true, emit: svLocusGraph_bin
    path "body.txt"                                                          , optional: true, emit: body_txt
    path "dragen.time_metrics.csv"                                           , optional: true, emit: dragen_time_metrics_csv
    path "header.txt"                                                        , optional: true, emit: header_txt
    path "match_log.small_indel_dedup.txt"                                   , optional: true, emit: match_log_small_indel_dedup_txt
    path "match_log.smn_dedup.txt"                                           , optional: true, emit: match_log_smn_dedup_txt
    path "streaming_log_*.csv"                                               , optional: true, emit: streaming_log_csv
    path "*_usage.txt"                                                       , emit: usage_txt
    path "**"                                                                , emit: all
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""

    // input data
    def bam_rex = /\.bam$/
    def cram_rex = /\.cram$/
    def fastq_rex = /\.f(ast)?q(\.gz)?$/
    def ora_rex = /\.ora$/

    // set RGID and RGSM - don't move this due to bug in version 24.10.4 and lower
    def fastq_ora = input =~ fastq_rex || input =~ ora_rex
    def rgid = args.contains("--RGID") || (!fastq_ora) ? "" : " --RGID ${meta.id}"
    def rgsm = args.contains("--RGSM") || (!fastq_ora) ? "" : " --RGSM ${meta.id}"
    args = args + rgid + rgsm


    def format_input = ""
    if (input instanceof List) {
        if (input.size() != 2) {
            error "Error: a maximum of 2 input files is supported."
        }
        if ((input[0] =~ fastq_rex && input[1] =~ fastq_rex) || (input[0] =~ ora_rex && input[1] =~ ora_rex)) {
            format_input = "-1 ${input[0]} -2 ${input[1]}"
        } else {
            error "Error: input of 2 files of this type is not supported."
        }
    } else {
        if (input =~ fastq_rex || input =~ ora_rex) {
            format_input = "-1 ${input}"
        } else if (input =~ bam_rex) {
            format_input = "-b ${input}"
        } else if (input =~ cram_rex) {
            format_input = "--cram-input ${input}"
        } else {
            error "Error: unsupported input type."
        }
    }

    // checkfingerprint-expected-vcf
    if (checkfingerprint_expected_vcf) {
        args = args + " --checkfingerprint-expected-vcf " + checkfingerprint_expected_vcf
    }

    // copy number variant calling
    if (cnv_combined_counts) {
        args = args + " --cnv-combined-counts " + cnv_combined_counts
    }
    if (cnv_exclude_bed) {
        args = args + " --cnv-exclude-bed " + cnv_exclude_bed
    }
    if (cnv_population_b_allele_vcf) {
        args = args + " --cnv-population-b-allele-vcf " + cnv_population_b_allele_vcf
    }
    if (cnv_segmentation_bed) {
        args = args + " --cnv-segmentation-bed " + cnv_segmentation_bed
    }
    if (cnv_target_bed) {
        args = args + " --cnv-target-bed " + cnv_target_bed
    }

    // cram-reference
    if (cram_reference) {
        args = args + " --cram-reference " + cram_reference
    }

    // dbsnp
    if (dbsnp) {
        args = args + " --dbsnp " + dbsnp
    }

    // fastqc
    if (fastqc_adapter_file) {
        args = args + " --fastqc-adapter-file " + fastqc_adapter_file
    }
    if (fastqc_kmer_file) {
        args = args + " --fastqc-kmer-file " + fastqc_kmer_file
    }

    // ora-reference
    if (ora_reference) {
        args = args + " --ora-reference " + ora_reference
    }

    // quality control
    if (qc_coverage_region) {
        if (qc_coverage_region.size() > 3) {
            error "Error: cannot have more than 3 qc-coverage-region files."
        } else {
            def region_args = qc_coverage_region.withIndex().collect { region, idx ->
                " --qc-coverage-region-${idx + 1} ${region}"
            }.join('')
            args = args + region_args
        }
    }
    if (qc_cross_cont_vcf) {
        args = args + " --qc-cross-cont-vcf " + qc_cross_cont_vcf
    }

    // ref-dir
    if (ref_dir) {
        args = args + " -r " + ref_dir
    }

    // repeat genotype
    if (repeat_genotype_ref_fasta) {
        args = args + " --repeat-genotype-ref-fasta " + repeat_genotype_ref_fasta
    }
    if (repeat_genotype_specs) {
        args = args + " --repeat-genotype-specs " + repeat_genotype_specs
    }

    // sample-sex
    if (meta.sex && !args.contains("--repeat-genotype-enable true")) {
        args = args + " --sample-sex " + meta.sex.toLowerCase()
    }

    // structural variation
    if (sv_call_regions_bed) {
        args = args + " --sv-call-regions-bed " + sv_call_regions_bed
    }
    if (sv_exclusion_bed) {
        args = args + " --sv-exclusion-bed " + sv_exclusion_bed
    }
    if (sv_forcegt_vcf) {
        args = args + " --sv-forcegt-vcf " + sv_forcegt_vcf
    }
    if (sv_systematic_noise) {
        args = args + " --sv-systematic-noise " + sv_systematic_noise
    }

    // trimming
    if (trim_adapter_read && trim_adapter_read.size() < 3) {
        args = args + " --trim-adapter-read1 " + trim_adapter_read.join(" --trim-adapter-read2 ")
    }
    if (trim_adapter_read_5prime && trim_adapter_read_5prime.size() < 3) {
        args = args + " --trim-adapter-r1-5prime " + trim_adapter_read_5prime.join(" --trim-adapter-r2-5prime ")
    }

    // variant-annotation-data
    if (variant_annotation_data) {
        args = args + " --variant-annotation-data " + variant_annotation_data
    }

    // variant calling
    if (vc_combine_phased_variants_distance_bed) {
        args = args + " --vc-combine-phased-variants-distance-bed " + vc_combine_phased_variants_distance_bed
    }
    if (vc_excluded_regions_bed) {
        args = args + " --vc-excluded-regions-bed " + vc_excluded_regions_bed
    }
    if (vc_forcegt_vcf) {
        args = args + " --vc-forcegt-vcf " + vc_forcegt_vcf
    }
    if (vc_log_bed) {
        args = args + " --vc-log-bed " + vc_log_bed
    }
    if (vc_mapping_metrics) {
        args = args + " --vc-mapping-metrics " + vc_mapping_metrics
    }
    if (vc_ml_dir) {
        args = args + " --vc-ml-dir " + vc_ml_dir
    }
    if (vc_ntd_error_params) {
        args = args + " --vc-ntd-error-params " + vc_ntd_error_params
    }
    if (vc_roh_blacklist_bed) {
        args = args + " --vc-roh-blacklist-bed " + vc_roh_blacklist_bed
    }
    if (vc_snp_error_cal_bed) {
        args = args + " --vc-snp-error-cal-bed " + vc_snp_error_cal_bed
    }
    if (vc_systematic_noise) {
        args = args + " --vc-systematic-noise " + vc_systematic_noise
    }
    if (vc_target_bed) {
        args = args + " --vc-target-bed " + vc_target_bed
    }

    // variant deduplication
    if (vd_eh_vcf) {
        args = args + " --vd-eh-vcf " + vd_eh_vcf
    }
    if (vd_small_variant_vcf) {
        args = args + " --vd-small-variant-vcf " + vd_small_variant_vcf
    }
    if (vd_sv_vcf) {
        args = args + " --vd-sv-vcf " + vd_sv_vcf
    }

    // vntr-catalog-bed
    if (vntr_catalog_bed) {
        args = args + " --vntr-catalog-bed " + vntr_catalog_bed
    }

    """
    dragen
        $args \\
        $format_input \\
        --output-file-prefix $prefix \\
        --output-directory \$(pwd)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$(dragen --version 2>&1) | sed 's/^dragen Version //;s/ Hash.*//')
    END_VERSIONS
    """

    stub:
    def VERSION = "stub"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch 20251120_usage.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: ${VERSION}
    END_VERSIONS
    """
}
