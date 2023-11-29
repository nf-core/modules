process DRAGEN {
    tag "$meta.id"
    label 'process_long'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'registry.hub.docker.com/etycksen/dragen4:4.2.4' }" // somehow the PATH is not containing dragen, when I run it locally it is there

    input:
    tuple val(meta), path(fastq), path(bam)
    path reference

    output:
    tuple val(meta), path("${prefix}.bam")                          , optional: true, emit: bam
    tuple val(meta), path("${prefix}.bai")                          , optional: true, emit: bai
    tuple val(meta), path("${prefix}.bam.md5sum")                   , optional: true, emit: bam_md5sum
    tuple val(meta), path("${prefix}.cnv.excluded_intervals.bed.gz"), optional: true, emit: cnv_excluded_intervals_bed_gz
    tuple val(meta), path("${prefix}.cnv.gff3")                     , optional: true, emit: cnv_gff3
    tuple val(meta), path("${prefix}.cnv.igv_session.xml")          , optional: true, emit: cnv_igv_session_xml
    tuple val(meta), path("${prefix}.cnv_metrics.csv")              , optional: true, emit: cnv_metrics_csv
    tuple val(meta), path("${prefix}.cnv.vcf")                      , optional: true, emit: cnv_vcf
    tuple val(meta), path("${prefix}.cnv.vcf.gz.md5sum")            , optional: true, emit: cnv_vcf_gz_md5sum
    tuple val(meta), path("${prefix}.cnv.vcf.gz.tbi")               , optional: true, emit: cnv_vcf_gz_tbi
    tuple val(meta), path("${prefix}.fastqc_metrics.csv")           , optional: true, emit: fastqc_metrics_csv
    tuple val(meta), path("${prefix}.fragment_length_hist.csv")     , optional: true, emit: fragment_length_hist_csv
    tuple val(meta), path("${prefix}.gc_metrics.csv")               , optional: true, emit: gc_metrics_csv
    tuple val(meta), path("${prefix}.g.vcf.gz")                     , optional: true, emit: g_vcf_gz
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi")                 , optional: true, emit: g_vcf_gz_tbi
    tuple val(meta), path("${prefix}.hard-filtered.gvcf.gz.md5sum") , optional: true, emit: hard_filtered_gvcf_gz_md5sum
    tuple val(meta), path("${prefix}.hard-filtered.baf.bw")         , optional: true, emit: hard_filtered_baf_bw
    tuple val(meta), path("${prefix}.hard-filtered.vcf")            , optional: true, emit: hard_filtered_vcf
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz")         , optional: true, emit: hard_filtered_vcf_gz
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz.md5sum")  , optional: true, emit: hard_filtered_vcf_gz_md5sum
    tuple val(meta), path("${prefix}.hard-filtered.vcf.gz.tbi")     , optional: true, emit: hard_filtered_vcf_gz_tbi
    tuple val(meta), path("${prefix}.improper.pairs.bw")            , optional: true, emit: improper_pairs_bw
    tuple val(meta), path("${prefix}.insert-stats.tab")             , optional: true, emit: insert_stats_tab
    tuple val(meta), path("${prefix}.mapping_metrics.csv")          , optional: true, emit: mapping_metrics_csv
    tuple val(meta), path("${prefix}.pcr-mode*.log")                , optional: true, emit: pcr_model_log
    tuple val(meta), path("${prefix}.ploidy.vcf")                   , optional: true, emit: ploidy_vcf
    tuple val(meta), path("${prefix}.ploidy.vcf.gz.md5sum")         , optional: true, emit: ploidy_vcf_gz_md5sum
    tuple val(meta), path("${prefix}.ploidy.vcf.gz.tbi")            , optional: true, emit: ploidy_vcf_gz_tbi
    tuple val(meta), path("${prefix}.ploidy_estimation_metrics.csv"), optional: true, emit: ploidy_estimation_metrics_csv
    tuple val(meta), path("${prefix}-replay.json")                  , optional: true, emit: replay_json
    tuple val(meta), path("${prefix}.repeats.bam")                  , optional: true, emit: repeats_bam
    tuple val(meta), path("${prefix}.repeats.vcf.gz")               , optional: true, emit: repeats_vcf_gz
    tuple val(meta), path("${prefix}.repeats.vcf.gz.tbi")           , optional: true, emit: repeats_vcf_gz_tbi
    tuple val(meta), path("${prefix}.roh.bed")                      , optional: true, emit: roh_bed
    tuple val(meta), path("${prefix}.roh_metrics.csv")              , optional: true, emit: roh_metrics_csv
    tuple val(meta), path("${prefix}.seg.bw")                       , optional: true, emit: seg_bw
    tuple val(meta), path("${prefix}.seg")                          , optional: true, emit: seg
    tuple val(meta), path("${prefix}.seg.c.partial.Y")              , optional: true, emit: seg_c_partial_y
    tuple val(meta), path("${prefix}.seg.called")                   , optional: true, emit: seg_called
    tuple val(meta), path("${prefix}.seg.called.merged")            , optional: true, emit: seg_called_merged
    tuple val(meta), path("${prefix}.sv.vcf")                       , optional: true, emit: sv_vcf
    tuple val(meta), path("${prefix}.sv.vcf.gz.tbi")                , optional: true, emit: sv_vcf_gz_tbi
    tuple val(meta), path("${prefix}.sv_metrics.csv")               , optional: true, emit: sv_metrics_csv
    tuple val(meta), path("${prefix}.target.counts.bw")             , optional: true, emit: target_counts_bw
    tuple val(meta), path("${prefix}.target.counts.diploid.bw")     , optional: true, emit: target_counts_diploid_bw
    tuple val(meta), path("${prefix}.target.counts.gz")             , optional: true, emit: target_counts_gz
    tuple val(meta), path("${prefix}.target.counts.gc-corrected.gz"), optional: true, emit: target_counts_gc_corrected_gz
    tuple val(meta), path("${prefix}.time_metrics.csv")             , optional: true, emit: time_metrics_csv
    tuple val(meta), path("${prefix}.trimmer_metrics.csv")          , optional: true, emit: trimmer_metrics_csv
    tuple val(meta), path("${prefix}.tn.bw")                        , optional: true, emit: tn_bw
    tuple val(meta), path("${prefix}.tn.tsv.gz")                    , optional: true, emit: tn_tsv_gz
    tuple val(meta), path("${prefix}.vc_hethom_ratio_metrics.csv")  , optional: true, emit: vc_hethom_ratio_metrics_csv
    tuple val(meta), path("${prefix}.vc_metrics.csv")               , optional: true, emit: vc_metrics_csv
    tuple val(meta), path("${prefix}.vcf.gz")                       , optional: true, emit: vcf_gz
    tuple val(meta), path("${prefix}.vcf.gz.md5sum")                , optional: true, emit: vcf_gz_md5sum
    tuple val(meta), path("${prefix}.vcf.gz.tbi")                   , optional: true, emit: vcf_gz_tbi
    tuple val(meta), path("${prefix}.wgs_fine_hist.csv")            , optional: true, emit: wgs_fine_hist_csv
    tuple val(meta), path("${prefix}.wgs_contig_mean_cov.csv")      , optional: true, emit: wgs_contig_mean_cov_csv
    tuple val(meta), path("${prefix}.wgs_coverage_metrics.csv")     , optional: true, emit: wgs_coverage_metrics_csv
    tuple val(meta), path("${prefix}.wgs_hist.csv")                 , optional: true, emit: wgs_hist_csv
    tuple val(meta), path("${prefix}.wgs_overall_mean_cov.csv")     , optional: true, emit: wgs_overall_mean_cov_csv
    path "./sort_spill/partitions.txt"                              , optional: true, emit: partitions_txt
    path "./sv/results/stats/alignmentStatsSummary.txt"             , optional: true, emit: alignmentStatsSummary_txt
    path "./sv/results/stats/candidate_metrics.csv"                 , optional: true, emit: candidate_metrics_csv
    path "./sv/results/stats/diploidSV.sv_metrics.csv"              , optional: true, emit: diploidSV_sv_metrics_csv
    path "./sv/results/stats/graph_metrics.csv"                     , optional: true, emit: graph_metrics_csv
    path "./sv/results/stats/svCandidateGenerationStats.tsv"        , optional: true, emit: svCandidateGenerationStats_tsv
    path "./sv/results/stats/svCandidateGenerationStats.xml"        , optional: true, emit: svCandidateGenerationStats_xml
    path "./sv/results/stats/svLocusGraphStats.tsv"                 , optional: true, emit: svLocusGraphStats_tsv
    path "./sv/results/variants/candidateSV.vcf.gz"                 , optional: true, emit: candidateSV_vcf_gz
    path "./sv/results/variants/candidateSV.vcf.gz.tbi"             , optional: true, emit: candidateSV_vcf_gz_tbi
    path "./sv/results/variants/diploidSV.vcf.gz"                   , optional: true, emit: diploidSV_vcf_gz
    path "./sv/results/variants/diploidSV.vcf.gz.tbi"               , optional: true, emit: diploidSV_vcf_gz_tbi
    path "./sv/workspace/alignmentStats.xml"                        , optional: true, emit: alignmentStats_xml
    path "./sv/workspace/chromDepth.txt"                            , optional: true, emit: chromDepth_txt
    path "./sv/workspace/edgeRuntimeLog.txt"                        , optional: true, emit: edgeRuntimeLog_txt
    path "./sv/workspace/genomeSegmentScanDebugInfo.txt"            , optional: true, emit: genomeSegmentScanDebugInfo_txt
    path "./sv/workspace/logs/config_log.txt"                       , optional: true, emit: config_log_txt
    path "./sv/workspace/svLocusGraph.bin"                          , optional: true, emit: svLocusGraph_bin
    path "body.txt"                                                 , emit: body_txt
    path "dragen.time_metrics.csv"                                  , emit: dragen_time_metrics_csv
    path "header.txt"                                               , emit: header_txt
    path "streaming_log_*.csv"                                      , emit: streaming_log_csv
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}" // change scope to global inside process
    def args = task.ext.args ?: ''
    def bin_path = task.ext.bin_path ?: 'dragen'

    def input = ''
    // from fastq or bam
    if (fastq) {
        if (fastq.size() > 2) {
            error "Error: cannot have more than 2 fastq files as input."
        } else {
            input = '-1 ' + fastq.join(' -2 ')
        }
    } else if (bam) {
        if (bam.size() > 1) {
            error "Error: cannot have more than 1 bam as input."
        } else {
            input = '-b ' + bam
        }
    }
    // set RGID and RGSM
    def rgid = args.contains("--RGID") ? "" : "--RGID ${meta.id}"
    def rgsm = args.contains("--RGSM") ? "" : "--RGSM ${meta.id}"

    """
    ${bin_path}_reset

    $bin_path \\
        $args \\
        $input \\
        $rgid \\
        $rgsm \\
        -n $task.cpus \\
        -r $reference \\
        --output-file-prefix $prefix \\
        --output-directory \$(pwd)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$($bin_path --version 2>&1) | sed 's/^dragen Version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$($bin_path --version 2>&1) | sed 's/^dragen Version //')
    END_VERSIONS
    """
}
