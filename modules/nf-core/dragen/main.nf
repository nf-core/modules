process DRAGEN {
    tag "$meta.id"
    label 'process_long'

    // no conda env or container image as Dragen requires specialized hardware to run

    // TODO: Problematic inputs
    //
    // build-sys-noise-vcfs-list: Text file containing the paths of normal VCFs. Specify the full VCF file paths. List one file per line.
    // cnv-normals-list: Specify text file that contains paths to the list of reference target counts files to be used as a panel of normals (new line separated).
    // config-file: Configuration file - could potentially set everything, no control
    // explify-sample-list: Path to .tsv file with sample names, associated FASTQs, etc.
    // SampleID		BatchID		RunID		ControlFlag		FastQs
    // MySample		MyBatch		MyRun		POS		/path/to/fastq1.gz		/path/to/fastq1.gz
    // fastq-list: Specifies CSV file that contains a list of FASTQ files to process.
    // fastq-list-sample-id: Specifies the sample ID for the list of FASTQ files specified by fastq-list.
    // ht-graph-vcf-list: Path to the text file containing the list of VCF files to be used to build the custom multigenome hash table
    // imputation-phase-input-list: Alternative to imputation-phase-input
    // imputation-phase-sample-type-list: Input file list of sample types where each line contains sample name NAME followed by TYPE. Required when imputing regions with ploidy that depends on this sample type.
    // input-batch-list: The path to a file containing a list of msVCF files to be merged, with the path to each file on a separate line. All the files listed must have been generated from the same global census file and all batches pertaining to that global census must be included in the merge.
    // input-census-list: File specifying list of census files for input (applicable when aggregate-censuses is true)
    // intermediate-results-dir: Specifies directory to store intermediate results in (eg, sort partitions).
    // ph-concat-all-input-list: PH concat all input msVCF file list (output files from phase rare step, in ascending position order), this option is exclusive with --ph-concat-all-input-list-sites-only
    // ph-concat-all-input-list-sites-only: Provides a .txt file with list of VCF containing all the haplotyped sites. The VCF files provided are the output files of Phase Rare step, in ascending position order, sex chromosomes at the end.
    // ph-ligate-common-input-list: PH ligate common input msVCF file list (output files from phase common step, in ascending position order)
    // ph-phase-common-input-list: PH phase common input msVCF file list
    // tumor-fastq-list: Inputs a CSV file containing a list of FASTQ files for the mapper, aligner, and somatic variant caller.
    // tumor-fastq-list-sample-id: Specifies the sample ID for the list of FASTQ files specified by tumor-fastq-list.
    //
    // TODO: meta.sex: male|female|none|auto ? found confounding flags --repeat-genotype-sex and --sample-sex; also it seems like it is used both with upper and lower case letters...

    input:
    tuple val(meta), path(fastq), path(ora_input) path(bam), path(cram), path(cnv_input), path(variant), path(dn_input_vcf), path(dn_sv_vcf), path(dn_cnv_vcf), path(imputation_phase_input), path(kmer_classifier_input_read_file), path(kmer_fq), path(maf_input_vcf), path(maf_input_json)
    tuple val(meta2), path(tumor_fastq), path(tumor_bam), path(tumor_cram)
    tuple val(meta3), path(bcl), path(samplesheet), path(run_info)
    tuple val(meta4), path(pedigree_cnv_input), path(pedigree_variant), path(pedigree_file)
    tuple val(meta5), path(input_census_file), path(input_cohort_file), path(input_global_census_file)
    tuple val(meta6), path(kmer_input_dir)
    tuple val(meta7), path(ph_phase_common_sample_type), path(ph_phase_common_scaffold), path(ph_phase_rare_input), path(ph_phase_rare_sample_type), path(ph_phase_rare_scaffold)
    tuple val(meta8), path(star_allele_cnv_vcf), path(star_allele_gvcf)
    path(amplicon_target_bed)
    path(annotation_file)
    path(atac_jaspar_database)
    path(checkfingerprint_expected_vcf)
    tuple path(cnv_combined_counts), path(cnv_exclude_bed), path(cnv_normal_b_allele_vcf), path(cnv_normal_cnv_vcf), path(cnv_normals_file), path(cnv_population_b_allele_vcf), path(cnv_segmentation_bed), path(cnv_somatic_essential_genes_bed), path(cnv_target_bed), path(cnv_within_gene_lr_bed)
    path(cram_reference)
    path(dbsnp)
    path(explify_ref_db_dir)
    tuple path(fastqc_adapter_file), path(fastqc_kmer_file)
    path(fragmentomics_exclude_bed)
    path(fragmentomics_wps_target_file)
    tuple path(gg_allele_list), path(gg_concurrency_regions), path(gg_output_regions), path(gg_regions), path(gg_regions_bed), path(gg_sites_list)
    tuple path(hrd_input_ascn), path(hrd_input_tn)
    tuple path(ht_alt_liftover), path(ht_decoys), path(ht_graph_bed), path(ht_hla_reference), path(ht_hla_reference), path(ht_mask_bed), path(ht_pop_alt_contigs), path(ht_pop_alt_liftover), path(ht_pop_snps), path(ht_reference)
    tuple path(imputation_chunk_input_region_list), path(imputation_genome_map_dir), path(imputation_ref_panel_dir)
    tuple path(kmer_classifier_db_file), path(kmer_classifier_db_to_taxid_json), path(kmer_classifier_reads_to_keep), path(kmer_detect_config_file), path(kmer_kmer_fasta), path(kmer_kmer_list)
    path(ml_recalibration_input_vcf)
    tuple path(msi_microsatellites_file), path(msi_ref_normal_dir)
    path(ora_reference)
    tuple path(ph_phase_common_config), path(ph_phase_common_map), path(ph_phase_common_reference), path(ph_phase_qc_estimation), path(ph_phase_qc_validation), path(ph_phase_rare_config), path(ph_phase_rare_map)
    tuple path(qc_coverage_region), path(qc_cross_cont_vcf), path(qc_somatic_contam_normal_pileup), path(qc_somatic_contam_tumor_pileup), path(qc_somatic_contam_vcf)
    path(ref_dir)
    tuple path(repeat_genotype_ref_fasta), path(repeat_genotype_specs)
    tuple path(rna_gf_blast_pairs), path(rna_gf_enriched_genes), path(rna_gf_enriched_regions), path(rna_gf_input_file), path(rna_repeat_genes), path(rna_repeat_intervals)
    tuple path(scatac_barcode_sequence_list), path(scatac_barcode_sequence_whitelist)
    tuple path(scrna_barcode_sequence_list), path(scrna_barcode_sequence_whitelist), path(scrna_cell_hashing_reference), path(scrna_demux_sample_vcf), path(scrna_demux_reference_vcf), path(scrna_feature_barcode_reference)
    path(single_cell_barcode_sequence_whitelist)
    tuple path(spatial_analysis_barcode_file), path(spatial_analysis_barcode_whitelist)
    tuple path(sv_call_regions_bed), path(sv_exclusion_bed), path(sv_forcegt_vcf), path(sv_locus_node_target_file), path(sv_somatic_ins_tandup_hotspot_regions_bed), path(sv_systematic_noise)
    path(tmb_ch_bed)
    tuple path(trim_adapter_read1), path(trim_adapter_read2)
    tuple path(trim_adapter_r1_5prime), path(trim_adapter_r2_5prime)
    tuple path(umi_correction_table), path(umi_fastq), path(umi_metrics_interval), path(umi_metrics_interval_file), path(umi_nonrandom_whitelist)
    path(variant_annotation_data)
    tuple path(vc_combine_phased_variants_distance_bed), path(vc_excluded_regions_bed), path(vc_forcegt_vcf), path(vc_log_bed), path(vc_mapping_metrics), path(vc_ml_dir), path(vc_ntd_error_params), path(vc_roh_blacklist_bed), path(vc_snp_error_cal_bed), path(vc_somatic_hotspots), path(vc_systematic_noise), path(vc_target_bed)
    tuple path(vd_eh_vcf), path(vd_small_variant_vcf), path(vd_sv_vcf)

    output:
    tuple val(meta), path("${prefix}.baf.bw") , optional: true, emit: .baf.bw
    tuple val(meta), path("${prefix}.bam")                                   , optional: true, emit: bam
    tuple val(meta), path("${prefix}.bai")                                   , optional: true, emit: bai
    tuple val(meta), path("${prefix}.bam.md5sum")                            , optional: true, emit: bam_md5sum
    tuple val(meta), path("${prefix}_chr_start-end.impute.phase.vcf.gz")     , optional: true, emit: chr_start_end_impute_phase_vcf_gz
    tuple val(meta), path("${prefix}.cnv.excluded_intervals.bed.gz")         , optional: true, emit: cnv_excluded_intervals_bed_gz
    tuple val(meta), path("${prefix}.cnv.gff3")                              , optional: true, emit: cnv_gff3
    tuple val(meta), path("${prefix}.cnv.igv_session.xml")                   , optional: true, emit: cnv_igv_session_xml
    tuple val(meta), path("${prefix}.cnv_metrics.csv")                       , optional: true, emit: cnv_metrics_csv
    tuple val(meta), path("${prefix}.cnv.vcf")                               , optional: true, emit: cnv_vcf
    tuple val(meta), path("${prefix}.cnv.vcf.gz")                            , optional: true, emit: cnv_vcf_gz
    tuple val(meta), path("${prefix}.cnv.vcf.gz.md5sum")                     , optional: true, emit: cnv_vcf_gz_md5sum
    tuple val(meta), path("${prefix}.cnv.vcf.gz.tbi")                        , optional: true, emit: cnv_vcf_gz_tbi
    tuple val(meta), path("${prefix}.fastqc_metrics.csv")                    , optional: true, emit: fastqc_metrics_csv
    tuple val(meta), path("${prefix}.fragment_length_hist.csv")              , optional: true, emit: fragment_length_hist_csv
    tuple val(meta), path("${prefix}.gc_metrics.csv")                        , optional: true, emit: gc_metrics_csv
    tuple val(meta), path("${prefix}.g.vcf.gz")                              , optional: true, emit: g_vcf_gz
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi")                          , optional: true, emit: g_vcf_gz_tbi
    tuple val(meta), path("${prefix}.hard-filtered.gvcf.gz.md5sum")          , optional: true, emit: hard_filtered_gvcf_gz_md5sum
    tuple val(meta), path("${prefix}.hard-filtered.baf.bw")                  , optional: true, emit: hard_filtered_baf_bw
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
    tuple val(meta), path("${prefix}.multiomics.barcodes.tsv.gz")            , optional: true, emit: multiomics_barcodes_tsv_gz
    tuple val(meta), path("${prefix}.multiomics.barcodeSummary.tsv")         , optional: true, emit: multiomics_barcodeSummary_tsv
    tuple val(meta), path("${prefix}.multiomics.features.tsv.gz")            , optional: true, emit: multiomics_features_tsv_gz
    tuple val(meta), path("${prefix}.multiomics.filtered.barcodes.tsv.gz")   , optional: true, emit: multiomics_filtered_barcodes.tsv_gz
    tuple val(meta), path("${prefix}.multiomics.matrix.mtx.gz")              , optional: true, emit: multiomics_matrix_mtx_gz
    tuple val(meta), path("${prefix}.multiomics.metrics.csv")                , optional: true, emit: multiomics_metrics_csv
    tuple val(meta), path("${prefix}.pcr-mode*.log")                         , optional: true, emit: pcr_model_log
    tuple val(meta), path("${prefix}.ploidy.vcf")                            , optional: true, emit: ploidy_vcf
    tuple val(meta), path("${prefix}.ploidy.vcf.gz.md5sum")                  , optional: true, emit: ploidy_vcf_gz_md5sum
    tuple val(meta), path("${prefix}.ploidy.vcf.gz.tbi")                     , optional: true, emit: ploidy_vcf_gz_tbi
    tuple val(meta), path("${prefix}.ploidy_estimation_metrics.csv")         , optional: true, emit: ploidy_estimation_metrics_csv
    tuple val(meta), path("${prefix}.preprocess.vcf.gz") , optional: true, emit: .preprocess.vcf.gz
    tuple val(meta), path("${prefix}.quant.genes.sf")                        , optional: true, emit: quant_genes_sf
    tuple val(meta), path("${prefix}.quant.metrics.csv")                     , optional: true, emit: quant_metrics_csv
    tuple val(meta), path("${prefix}.quant.sf")                              , optional: true, emit: quant_sf
    tuple val(meta), path("${prefix}.quant.transcript_coverage.txt")         , optional: true, emit: quant_transcript_coverage_txt
    tuple val(meta), path("${prefix}.quant.transcript_fragment_lengths.txt") , optional: true, emit: quant_transcript_fragment_lengths_txt
    tuple val(meta), path("${prefix}-replay.json")                           , optional: true, emit: replay_json
    tuple val(meta), path("${prefix}.repeats.bam")                           , optional: true, emit: repeats_bam
    tuple val(meta), path("${prefix}.repeats.vcf.gz")                        , optional: true, emit: repeats_vcf_gz
    tuple val(meta), path("${prefix}.repeats.vcf.gz.tbi")                    , optional: true, emit: repeats_vcf_gz_tbi
    tuple val(meta), path("${prefix}.roh.bed")                               , optional: true, emit: roh_bed
    tuple val(meta), path("${prefix}.roh_metrics.csv")                       , optional: true, emit: roh_metrics_csv
    tuple val(meta), path("${prefix}.scATAC.barcodeSummary.tsv")             , optional: true, emit: scATAC_barcodeSummary_tsv
    tuple val(meta), path("${prefix}.scATAC.barcodes.tsv")                   , optional: true, emit: scATAC_barcodes_tsv
    tuple val(meta), path("${prefix}.scATAC.barcodes.tsv.gv")                , optional: true, emit: scATAC_barcodes_tsv_gv
    tuple val(meta), path("${prefix}.scATAC.filtered.barcodes.tsv.gz")       , optional: true, emit: scATAC_filtered_barcodes_tsv_gz
    tuple val(meta), path("${prefix}.scATAC.filtered.matrix.mtx.gz")         , optional: true, emit: scATAC_filtered_matrix_mtx_gz
    tuple val(meta), path("${prefix}.scATAC.fragments.bigwig")               , optional: true, emit: scATAC_fragments_bigwig
    tuple val(meta), path("${prefix}.scATAC.fragments.tsv")                  , optional: true, emit: scATAC_fragments_tsv
    tuple val(meta), path("${prefix}.scATAC.matrix.mtx.gv")                  , optional: true, emit: scATAC_matrix_mtx_gv
    tuple val(meta), path("${prefix}.scATAC.metrics.csv")                    , optional: true, emit: scATAC_metrics_csv
    tuple val(meta), path("${prefix}.scATAC.motifs.matrix.mtx")              , optional: true, emit: scATAC_motifs_matrix_mtx
    tuple val(meta), path("${prefix}.scATAC.peaks.tsv")                      , optional: true, emit: scATAC_peaks_tsv
    tuple val(meta), path("${prefix}.scATAC.peaks.tsv.gv")                   , optional: true, emit: scATAC_peaks_tsv_gv
    tuple val(meta), path("${prefix}.scATAC.tf.motifs.tsv")                  , optional: true, emit: scATAC_tf_motifs_tsv
    tuple val(meta), path("${prefix}.scRNA.barcodeSummary.tsv")              , optional: true, emit: scRNA_barcodeSummary_tsv
    tuple val(meta), path("${prefix}.scRNA.barcodes.tsv.gz")                 , optional: true, emit: scRNA_barcodes_tsv_gz
    tuple val(meta), path("${prefix}.scRNA.demux.tsv")                       , optional: true, emit: scRNA_demux_tsv
    tuple val(meta), path("${prefix}.scRNA.filtered.barcodes.tsv.gz")        , optional: true, emit: scRNA_filtered_barcodes_tsv_gz
    tuple val(meta), path("${prefix}.scRNA.filtered.matrix.mtx.gz")          , optional: true, emit: scRNA_filtered_matrix_mtx_gz
    tuple val(meta), path("${prefix}.scRNA.genes.tsv.gz")                    , optional: true, emit: scRNA_genes_tsv_gz
    tuple val(meta), path("${prefix}.scRNA.matrix.mtx.gz")                   , optional: true, emit: scRNA_matrix_mtx_gz
    tuple val(meta), path("${prefix}.scRNA.metrics.csv")                     , optional: true, emit: scRNA_metrics_csv
    tuple val(meta), path("${prefix}.scRNA.metrics.demuxSamples.csv")        , optional: true, emit: scRNA_metrics_demuxSamples_csv
    tuple val(meta), path("${prefix}.seg.bw")                                , optional: true, emit: seg_bw
    tuple val(meta), path("${prefix}.seg")                                   , optional: true, emit: seg
    tuple val(meta), path("${prefix}.seg.c.partial.Y")                       , optional: true, emit: seg_c_partial_y
    tuple val(meta), path("${prefix}.seg.called")                            , optional: true, emit: seg_called
    tuple val(meta), path("${prefix}.seg.called.merged")                     , optional: true, emit: seg_called_merged
    tuple val(meta), path("${prefix}.select.gvcf")                           , optional: true, emit: select_gvcf
    tuple val(meta), path("${prefix}.SJ.saturation.txt")                     , optional: true, emit: SJ_saturation_txt
    tuple val(meta), path("${prefix}.small_indel_dedup")                     , optional: true, emit: small_indel_dedup
    tuple val(meta), path("${prefix}.smn_dedup")                             , optional: true, emit: smn_dedup
    tuple val(meta), path("${prefix}.star_allele.tsv")                       , optional: true, emit: star_allele_tsv
    tuple val(meta), path("${prefix}.star_allele.json")                      , optional: true, emit: star_allele_json
    tuple val(meta), path("${prefix}.sv.vcf")                                , optional: true, emit: sv_vcf
    tuple val(meta), path("${prefix}.sv.vcf.gz.tbi")                         , optional: true, emit: sv_vcf_gz_tbi
    tuple val(meta), path("${prefix}.sv_metrics.csv")                        , optional: true, emit: sv_metrics_csv
    tuple val(meta), path("${prefix}.target.counts.bw")                      , optional: true, emit: target_counts_bw
    tuple val(meta), path("${prefix}.target.counts.diploid.bw")              , optional: true, emit: target_counts_diploid_bw
    tuple val(meta), path("${prefix}.target.counts.gz")                      , optional: true, emit: target_counts_gz
    tuple val(meta), path("${prefix}.target.counts.gc-corrected.gz")         , optional: true, emit: target_counts_gc_corrected_gz
    tuple val(meta), path("${prefix}.targeted.json")                         , optional: true, emit: targeted_json
    tuple val(meta), path("${prefix}.targeted.vcf.gz")                       , optional: true, emit: targeted_vcf_gz
    tuple val(meta), path("${prefix}.time_metrics.csv")                      , optional: true, emit: time_metrics_csv
    tuple val(meta), path("${prefix}.trimmer_metrics.csv")                   , optional: true, emit: trimmer_metrics_csv
    tuple val(meta), path("${prefix}.tn.bw")                                 , optional: true, emit: tn_bw
    tuple val(meta), path("${prefix}.tn.tsv.gz")                             , optional: true, emit: tn_tsv_gz
    tuple val(meta), path("${prefix}.vc_hethom_ratio_metrics.csv")           , optional: true, emit: vc_hethom_ratio_metrics_csv
    tuple val(meta), path("${prefix}.vc_metrics.csv")                        , optional: true, emit: vc_metrics_csv
    tuple val(meta), path("${prefix}.vcf.gz")                                , optional: true, emit: vcf_gz
    tuple val(meta), path("${prefix}.vcf.gz.md5sum")                         , optional: true, emit: vcf_gz_md5sum
    tuple val(meta), path("${prefix}.vcf.gz.tbi")                            , optional: true, emit: vcf_gz_tbi
    tuple val(meta), path("${prefix}.wgs_fine_hist.csv")                     , optional: true, emit: wgs_fine_hist_csv
    tuple val(meta), path("${prefix}.wgs_contig_mean_cov.csv")               , optional: true, emit: wgs_contig_mean_cov_csv
    tuple val(meta), path("${prefix}.wgs_coverage_metrics.csv")              , optional: true, emit: wgs_coverage_metrics_csv
    tuple val(meta), path("${prefix}.wgs_hist.csv")                          , optional: true, emit: wgs_hist_csv
    tuple val(meta), path("${prefix}.wgs_overall_mean_cov.csv")              , optional: true, emit: wgs_overall_mean_cov_csv
    path "./sort_spill/partitions.txt"                                       , optional: true, emit: partitions_txt
    path "./sv/results/stats/alignmentStatsSummary.txt"                      , optional: true, emit: alignmentStatsSummary_txt
    path "./sv/results/stats/candidate_metrics.csv"                          , optional: true, emit: candidate_metrics_csv
    path "./sv/results/stats/diploidSV.sv_metrics.csv"                       , optional: true, emit: diploidSV_sv_metrics_csv
    path "./sv/results/stats/graph_metrics.csv"                              , optional: true, emit: graph_metrics_csv
    path "./sv/results/stats/somaticSV.sv_metrics.csv"                       , optional: true, emit: somaticSV_sv_metrics_csv
    path "./sv/results/stats/svCandidateGenerationStats.tsv"                 , optional: true, emit: svCandidateGenerationStats_tsv
    path "./sv/results/stats/svCandidateGenerationStats.xml"                 , optional: true, emit: svCandidateGenerationStats_xml
    path "./sv/results/stats/svLocusGraphStats.tsv"                          , optional: true, emit: svLocusGraphStats_tsv
    path "./sv/results/stats/tumorSV.sv_metrics.csv"                         , optional: true, emit: tumorSV_sv_metrics_csv
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
    path "match_log.small_indel_dedup.txt"                                   , optional: true, mit: match_log_small_indel_dedup_txt
    path "match_log.smn_dedup.txt"                                           , optional: true, mit: match_log_smn_dedup_txt
    path "*_usage.txt"                                                       , emit: usage_txt
    path "body.txt"                                                          , emit: body_txt
    path "dragen.time_metrics.csv"                                           , emit: dragen_time_metrics_csv
    path "header.txt"                                                        , emit: header_txt
    path "streaming_log_*.csv"                                               , emit: streaming_log_csv
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}" // change scope to global inside process
    def args = task.ext.args ?: ""

    // input data
    def input = ""
    if (bcl) {
        input = "--bcl-input-directory " + bcl
        if (samplesheet) {
            input = input + " --sample-sheet " + samplesheet
        }
        if (run_info) {
            input = input + " --run-info " + run_info
        }
    } else if (fastq || tumor_fastq) {
        if (fastq.size() > 2 || tumor_fastq.size() > 2) {
            error "Error: cannot have more than 2 fastq files as input."
        } else {
            if (tumor_fastq) {
                input = "--tumor-fastq1 " + tumor_fastq.join(" --tumor-fastq2 ")
            }
            if (fastq) {
                input = input + " -1 " + fastq.join(" -2 ")
            }
        }
    } else if (ora_input) {
        if (ora_input.size() > 2) {
            error "Error: cannot have more than 2 ora files as input."
        } else {
            input = "--ora-input " + ora_input.join(" --ora-input2 ")
        }
    } else if (bam || tumor_bam) {
        if (bam.size() > 1 || tumor_bam.size() > 1) {
            error "Error: cannot have more than 1 bam as input."
        } else {
            if (tumor_bam) {
                input = "--tumor-bam-input " + tumor_bam
            }
            if (bam) {
                input = input + " -b " + bam
            }
        }
    } else if (cram || tumor_cram) {
        if (cram.size() > 1 || tumor_cram.size() > 1) {
            error "Error: cannot have more than 1 cram as input."
        } else {
            if (tumor_cram) {
                input = "--tumor-cram-input " + tumor_cram
            }
            if (cram) {
                input = input + " --cram-input " + cram
            }
        }
    } else if (cnv_input) {
        input = "--cnv-input " + cnv_input
        if (pedigree_file) {
            input = input + " --pedigree-file " + pedigree_file + " --cnv-input " + pedigree_cnv_input.join(" --cnv-input ")
        }
    } else if (variant) {
        input = " --variant " + variant.join(" --variant ")
        if (pedigree_file) {
            input = input + " --pedigree-file " + pedigree_file + " --variant " + pedigree_variant.join(" --variant ")
        }
    } else if (dn_input_vcf) { // do cnv and sv only occur if input vcf is also provided or can they occur on their own
        input = " --dn-input-vcf " + dn_input_vcf + " --dn-output-vcf " + prefix + ".dn_output.vcf"
        if (dn_sv_vcf) {
            input = input + " --dn-sv-vcf " + dn_sv_vcf
        }
        if (dn_cnv_vcf) {
            input = input + " --dn-cnv-vcf " + dn_cnv_vcf
        }
    } else if (imputation_phase_input) {
        input = " --imputation-phase-input " + imputation_phase_input + " --imputation-phase-sample-type " + meta.imputation_type
    } else if (input_census_file && input_cohort_file && input_global_census_file) {
        input = " --input-census-file " + input_census_file + " --input-cohort-file " + input_cohort_file + " --input-global-census-file " + input_global_census_file
    } else if (kmer_classifier_input_read_file) {
        input = " --kmer-classifier-input-read-file " + kmer_classifier_input_read_file
    } else if (kmer_fq) {
        if (kmer_fq.size() > 2) {
            error "Error: cannot have more than 2 fastq as input."
        } else {
            input = " --kmer-fq " + kmer_fq.join(" --kmer-fq2 ")
        }
    } else if (kmer_input_dir) {
        input = " --kmer-input-dir " + kmer_input_dir
    } else if (maf_input_vcf && maf_input_json) { // both or only one?
        input = " --maf-input-vcf " + maf_input_vcf + " --maf-input-json " + maf_input_json
    } else if (ph_phase_common_sample_type) {
        input = " --ph-phase-common-sample-type " + ph_phase_common_sample_type
        if (ph_phase_common_scaffold) {
            input = input + " --ph-phase-common-scaffold " + ph_phase_common_scaffold
        }
    } else if (ph_phase_rare_input) {
        input = " --ph-phase-rare-input " + ph_phase_rare_input + " --ph-phase-rare-sample-type " + ph_phase_rare_sample_type +  " --ph-phase-rare-scaffold " + ph_phase_rare_scaffold
    } else if (star_allele_cnv_vcf && star_allele_gvcf) {
        input = " --star-allele-cnv-vcf " + star_allele_cnv_vcf + " --star-allele-gvcf " + star_allele_gvcf
    } else {
        error "Error: missing input."
    }

    // amplicon-target-bed
    if (amplicon_target_bed) {
        args = args + " --amplicon-target-bed " + amplicon_target_bed
    }

    // annotation-file
    if (annotation_file) {
        args = args + " -a " + annotation_file
    }

    // atac-jaspar-database
    if (atac_jaspar_database) {
        args = args + " --atac-jaspar-database " + atac_jaspar_database
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
    if (cnv_normal_b_allele_vcf) {
        args = args + " --cnv-normal-b-allele-vcf " + cnv_normal_b_allele_vcf
    }
    if (cnv_normal_cnv_vcf) {
        args = args + " --cnv-normal-cnv-vcf " + cnv_normal_cnv_vcf
    }
    if (cnv_normals_file) {
        args = args + " --cnv-normals-file " + cnv_normals_file.join(" --cnv-normals-file ")
    }
    if (cnv_population_b_allele_vcf) {
        args = args + " --cnv-population-b-allele-vcf " + cnv_population_b_allele_vcf
    }
    if (cnv_segmentation_bed) {
        args = args + " --cnv-segmentation-bed " + cnv_segmentation_bed
    }
    if (cnv_somatic_essential_genes_bed) {
        args = args + " --cnv-somatic-essential-genes-bed " + cnv_somatic_essential_genes_bed
    }
    if (cnv_target_bed) {
        args = args + " --cnv-target-bed " + cnv_target_bed
    }
    if (cnv_within_gene_lr_bed) {
        args = args + " --cnv-within-gene-lr-bed " + cnv_within_gene_lr_bed
    }

    // cram-reference
    if (cram_reference) {
        args = args + " --cram-reference " + cram_reference
    }

    // dbsnp
    if (dbsnp) {
        args = args + " --dbsnp " + dbsnp
    }

    // explify-ref-db-dir
    if (explify_ref_db_dir) {
        args = args + " --explify-ref-db-dir " + explify_ref_db_dir
    }

    // fastqc
    if (fastqc_adapter_file) {
        args = args + " --fastqc-adapter-file " + fastqc_adapter_file
    }
    if (fastqc_kmer_file) {
        args = args + " --fastqc-kmer-file " + fastqc_kmer_file
    }

    // fragmentomics
    if (fragmentomics_exclude_bed) {
        args = args + " --fragmentomics-exclude-bed " + fragmentomics_exclude_bed
    }
    if (fragmentomics_wps_target_file) {
        args = args + " --fragmentomics-wps-target-file " + fragmentomics_wps_target_file
    }

    // gvcf genotyper analyser
    if (gg_allele_list) {
        args = args + " --gg-allele-list " + gg_allele_list
    }
    if (gg_concurrency_regions) {
        args = args + " --gg-concurrency-regions " + gg_concurrency_regions
    }
    if (gg_output_regions) {
        args = args + " --gg-output-regions " + gg_output_regions
    }
    if (gg_regions) {
        args = args + " --gg-regions " + gg_regions
    }
    if (gg_regions_bed) {
        args = args + " --gg-regions-bed " + gg_regions_bed
    }
    if (gg_sites_list) {
        args = args + " --gg-sites-list " + gg_sites_list
    }

    // homologous recombination deficiency
    if (hrd_input_tn) {
        args = args + " --hrd-input-tn " + hrd_input_tn
    }
    if (hrd_input_ascn) {
        args = args + " --hrd-input-ascn " + hrd_input_ascn
    }

    // hash table
    if (ht_alt_liftover) {
        args = args + " --ht-alt-liftover " + ht_alt_liftover
    }
    if (ht_decoys) {
        args = args + " --ht-decoys " + ht_decoys
    }
    if (ht_graph_bed) {
        args = args + " --ht-graph-bed " + ht_graph_bed
    }
    if (ht_hla_reference) {
        args = args + " --ht-hla-reference " + ht_hla_reference
    }
    if (ht_hla_reference) {
        args = args + " --ht-hla-reference " + ht_hla_reference
    }
    if (ht_mask_bed) {
        args = args + " --ht-mask-bed " + ht_mask_bed
    }
    if (ht_pop_alt_contigs) {
        args = args + " --ht-pop-alt-contigs " + ht_pop_alt_contigs
    }
    if (ht_pop_alt_liftover) {
        args = args + " --ht-pop-alt-liftover " + ht_pop_alt_liftover
    }
    if (ht_pop_snps) {
        args = args + " --ht-pop-snps " + ht_pop_snps
    }
    if (ht_reference) {
        args = args + " --ht-reference " + ht_reference
    }

    // imputation
    if (imputation_chunk_input_region_list) {
        args = args + " --imputation-chunk-input-region-list " + imputation_chunk_input_region_list
    }
    if (imputation_genome_map_dir) {
        args = args + " --imputation-genome-map-dir " + imputation_genome_map_dir
    }
    if (imputation_ref_panel_dir) {
        args = args + " --imputation-ref-panel-dir " + imputation_ref_panel_dir
    }

    // kmer
    if (kmer_classifier_db_file) {
        args = args + " --kmer-classifier-db-file " + kmer_classifier_db_file
    }
    if (kmer_classifier_db_to_taxid_json) {
        args = args + " --kmer-classifier-db-to-taxid-json " + kmer_classifier_db_to_taxid_json
    }
    if (kmer_classifier_reads_to_keep) {
        args = args + " --kmer-classifier-reads-to-keep " + kmer_classifier_reads_to_keep
    }
    if (kmer_detect_config_file) {
        args = args + " --kmer-detect-config-file " + kmer_detect_config_file
    }
    if (kmer_kmer_fasta) {
        args = args + " --kmer-kmer-fasta " + kmer_kmer_fasta
    }
    if (kmer_kmer_list) {
        args = args + " --kmer-kmer-list " + kmer_kmer_list
    }

    // ml-recalibration-input-vcf
    if (ml_recalibration_input_vcf) {
        args = args + " --ml-recalibration-input-vcf " + ml_recalibration_input_vcf
    }

    // micro-satellite instability
    if (msi_microsatellites_file) {
        args = args + " --msi-microsatellites-file " + msi_microsatellites_file
    }
    if (msi_ref_normal_dir) {
        args = args + " --msi-ref-normal-dir " + msi_ref_normal_dir
    }

    // ora-reference
    if (ora_reference) {
        args = args + " --ora-reference " + ora_reference
    }

    // population haplotyping
    if (ph_phase_common_config) {
        args = args + " --ph-phase-common-config " + ph_phase_common_config
    }
    if (ph_phase_common_map) {
        args = args + " --ph-phase-common-map " + ph_phase_common_map
    }
    if (ph_phase_common_reference) {
        args = args + " --ph-phase-common-reference " + ph_phase_common_reference
    }
    if (ph_phase_qc_estimation) {
        args = args + " --ph-phase-qc-estimation " + ph_phase_qc_estimation
    }
    if (ph_phase_qc_validation) {
        args = args + " --ph-phase-qc-validation " + ph_phase_qc_validation
    }
    if (ph_phase_rare_config) {
        args = args + " --ph-phase-rare-config " + ph_phase_rare_config
    }
    if (ph_phase_rare_map) {
        args = args + " --ph-phase-rare-map " + ph_phase_rare_map
    }

    // quality control
    if (qc_coverage_region) {
        if (qc_coverage_region.size() > 3) {
            error "Error: cannot have more than 3 qc-coverage-region files."
        } else {
            for (int i = 1; i < qc_coverage_region.size() + 1; i++) {
                args = args + " --qc-coverage-region-" + i + " " + qc_coverage_region[i]
            }
        }
    }
    if (qc_cross_cont_vcf) {
        args = args + " --qc-cross-cont-vcf " + qc_cross_cont_vcf
    }
    if (qc_somatic_contam_normal_pileup) {
        args = args + " --qc-somatic-contam-normal-pileup " + qc_somatic_contam_normal_pileup
    }
    if (qc_somatic_contam_tumor_pileup) {
        args = args + " --qc-somatic-contam-tumor-pileup " + qc_somatic_contam_tumor_pileup
    }
    if (qc_somatic_contam_vcf) {
        args = args + " --qc-somatic-contam-vcf " + qc_somatic_contam_vcf
    }

    // ref-dir
    if (ref_dir) {
        args = args + " -r " + ref_dir
    }

    // repeat genotype
    if (args.contains("--repeat-genotype-enable true")) {
        // TODO: check if all required input is defined
        args = args + " --repeat-genotype-sex " + ${meta.sex}.toLowerCase() + " --repeat-genotype-ref-fasta " + repeat_genotype_ref_fasta + " --repeat-genotype-specs " + repeat_genotype_specs
    }

    // rna
    if (rna_gf_blast_pairs) {
        args = args + " --rna-gf-blast-pairs " + rna_gf_blast_pairs
    }
    if (rna_gf_enriched_genes) {
        args = args + " --rna-gf-enriched-genes " + rna_gf_enriched_genes
    }
    if (rna_gf_enriched_regions) {
        args = args + " --rna-gf-enriched-regions " + rna_gf_enriched_regions
    }
    if (rna_gf_input_file) {
        args = args + " --rna-gf-input-file " + rna_gf_input_file
    }
    if (rna_repeat_genes) {
        args = args + " --rna-repeat-genes " + rna_repeat_genes
    }
    if (rna_repeat_intervals) {
        args = args + " --rna-repeat-intervals " + rna_repeat_intervals
    }

    // single cell atac
    if (scatac_barcode_sequence_list) {
        args = args + " --scatac-barcode-sequence-list " + scatac_barcode_sequence_list
    }
    if (scatac_barcode_sequence_whitelist) {
        args = args + " --scatac-barcode-sequence-whitelist " + scatac_barcode_sequence_whitelist
    }

    // single cell rna
    if (scrna_barcode_sequence_list) {
        args = args + " --scrna-barcode-sequence-list " + scrna_barcode_sequence_list
    }
    if (scrna_barcode_sequence_whitelist) {
        args = args + " --scrna-barcode-sequence-whitelist " + scrna_barcode_sequence_whitelist
    }
    if (scrna_cell_hashing_reference) {
        args = args + " --scrna-cell-hashing-reference " + scrna_cell_hashing_reference
    }
    if (scrna_demux_sample_vcf) {
        args = args + " --scrna-demux-sample-vcf " + scrna_demux_sample_vcf
    }
    if (scrna_demux_reference_vcf) {
        args = args + " --scrna-demux-reference-vcf " + scrna_demux_reference_vcf
    }
    if (scrna_feature_barcode_reference) {
        args = args + " --scrna-feature-barcode-reference " + scrna_feature_barcode_reference
    }

    // single-cell-barcode-sequence-whitelist
    if (single_cell_barcode_sequence_whitelist) {
        args = args + " --single-cell-barcode-sequence-whitelist " + single_cell_barcode_sequence_whitelist
    }

    // spatial analysis
    if (spatial_analysis_barcode_file) {
        args = args + " --spatial-analysis-barcode-file " + spatial_analysis_barcode_file
    }
    if (spatial_analysis_barcode_whitelist) {
        args = args + " --spatial-analysis-barcode-whitelist " + spatial_analysis_barcode_whitelist
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
    if (sv_locus_node_target_file) {
        args = args + " --sv-locus-node-target-file " + sv_locus_node_target_file
    }
    if (sv_somatic_ins_tandup_hotspot_regions_bed) {
        args = args + " --sv-somatic-ins-tandup-hotspot-regions-bed " + sv_somatic_ins_tandup_hotspot_regions_bed
    }
    if (sv_systematic_noise) {
        args = args + " --sv-systematic-noise " + sv_systematic_noise
    }

    // tmb-ch-bed
    if (tmb_ch_bed) {
        args = args + " --tmb-ch-bed " + tmb_ch_bed
    }

    // trim-adapter-read
    if (trim_adapter_read1 && trim_adapter_read2) {
        args = args + " --trim-adapter-read1 " + trim_adapter_read1 + " --trim-adapter-read2 " + trim_adapter_read2
    }

    // trim-adapter-read-5prime
    if (trim_adapter_r1_5prime && trim_adapter_r2_5prime) {
        args = args + " --trim-adapter-r1-5prime " + trim_adapter_r1_5prime + " --trim-adapter-r2-5prime " + trim_adapter_r2_5prime
    }

    // unique molecular identifier
    if (umi_correction_table) {
        args = args + " --umi-correction-table " + umi_correction_table
    }
    if (umi_fastq) {
        args = args + " --umi-fastq " + umi_fastq
    }
    if (umi_metrics_interval) {
        args = args + " --umi-metrics-interval " + umi_metrics_interval
    }
    if (umi_metrics_interval_file) {
        args = args + " --umi-metrics-interval-file " + umi_metrics_interval_file
    }
    if (umi_nonrandom_whitelist) {
        args = args + " --umi-nonrandom-whitelist " + umi_nonrandom_whitelist
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
    if (vc_somatic_hotspots) {
        args = args + " --vc-somatic-hotspots " + vc_somatic_hotspots
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

    // set RGID and RGSM
    def rgid = args.contains("--RGID") || (!fastq && !ora_input) ? "" : " --RGID ${meta.id}"
    def rgid_tumor = args.contains("--RGID-tumor") || !tumor_fastq ? "" : " --RGID-tumor ${meta.id}"
    def rgsm = args.contains("--RGSM") || (!fastq && !ora_input) ? "" : " --RGSM ${meta.id}"
    def rgsm_tumor = args.contains("--RGSM-tumor") || !tumor_fastq ? "" : " --RGSM-tumor ${meta.id}"
    args = args + rgid + rgid_tumor + rgsm + rgsm_tumor

    """
    dragen
        $args \\
        $input \\
        -n $task.cpus \\
        --output-file-prefix $prefix \\
        --output-directory \$(pwd)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$(dragen --version 2>&1) | sed 's/^dragen Version //')
    END_VERSION
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$(dragen --version 2>&1) | sed 's/^dragen Version //')
    END_VERSIONS
    """
}
