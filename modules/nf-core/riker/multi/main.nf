process RIKER_MULTI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/54/54c7820b49cfb5fada32c1825ac8a46c05d0085105c50055cbce4c37701ab95e/data' :
        'community.wave.seqera.io/library/riker:0.4.0--4e7eeb0beed906c0' }"

    input:
    tuple val(meta),  path(bam), path(bai), path(error_vcf), path(error_vcf_idx), path(error_intervals), path(gcbias_exclude_intervals), path(hybcap_baits, stageAs: 'baits/*'), path(hybcap_targets, stageAs: 'targets/*'), path(rna_gene_model), path(rna_ribosomal_intervals), path(wgs_intervals)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.alignment-metrics.txt"),          emit: alignment_metrics,         optional: true
    tuple val(meta), path("*.base-distribution-by-cycle.txt"), emit: base_dist,                 optional: true
    tuple val(meta), path("*.error-indel.txt"),                emit: error_indel,               optional: true
    tuple val(meta), path("*.error-mismatch.txt"),             emit: error_mismatch,            optional: true
    tuple val(meta), path("*.error-overlap.txt"),              emit: error_overlap,             optional: true
    tuple val(meta), path("*.gcbias-detail.txt"),              emit: gcbias_detail,             optional: true
    tuple val(meta), path("*.gcbias-summary.txt"),             emit: gcbias_summary,            optional: true
    tuple val(meta), path("*.hybcap-metrics.txt"),             emit: hybcap_metrics,            optional: true
    tuple val(meta), path("*.hybcap-per-base.txt*"),           emit: hybcap_per_base,           optional: true
    tuple val(meta), path("*.hybcap-per-target.txt"),          emit: hybcap_per_target,         optional: true
    tuple val(meta), path("*.isize-histogram.txt"),            emit: isize_histogram,           optional: true
    tuple val(meta), path("*.isize-metrics.txt"),              emit: isize_metrics,             optional: true
    tuple val(meta), path("*.mean-quality-by-cycle.txt"),      emit: mean_qual,                 optional: true
    tuple val(meta), path("*.pdf"),                            emit: pdf,                       optional: true
    tuple val(meta), path("*.quality-score-distribution.txt"), emit: qual_dist,                 optional: true
    tuple val(meta), path("*.rna-biotype.txt"),                emit: rna_biotype,               optional: true
    tuple val(meta), path("*.rna-insert-size-histogram.txt"),  emit: rna_insert_size_histogram, optional: true
    tuple val(meta), path("*.rna-insert-size.txt"),            emit: rna_insert_size,           optional: true
    tuple val(meta), path("*.rna-metrics.txt"),                emit: rna_metrics,               optional: true
    tuple val(meta), path("*.wgs-coverage.txt"),               emit: wgs_coverage,              optional: true
    tuple val(meta), path("*.wgs-metrics.txt"),                emit: wgs_metrics,               optional: true
    tuple val("${task.process}"), val('riker'), eval("riker --version 2>&1 | sed 's/riker //'") , topic: versions, emit: versions_riker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def reference_arg = fasta ? "--reference ${fasta}" : ''
    if ((hybcap_baits as Boolean) ^ (hybcap_targets as Boolean)) {
        error "RIKER_MULTI: both 'baits' and 'targets' must be provided together, or neither"
    }
    def error_vcf_arg = error_vcf && error_vcf_idx ? "--error::vcf ${error_vcf}" : ''
    def error_intervals_arg = error_intervals ? "--error::intervals ${error_intervals}" : ''
    def gcbias_exclude_intervals_arg = gcbias_exclude_intervals ? "--gcbias::exclude-intervals ${gcbias_exclude_intervals}" : ''
    def hybcap_opts = (hybcap_baits && hybcap_targets) ? "--hybcap::baits ${hybcap_baits} --hybcap::targets ${hybcap_targets}" : ''
    def rna_gene_model_arg = rna_gene_model ? "--rna::gene-model ${rna_gene_model}" : ''
    def rna_ribosomal_intervals_arg = rna_ribosomal_intervals ? "--rna::ribosomal-intervals ${rna_ribosomal_intervals}" : ''
    def wgs_intervals_arg = wgs_intervals ? "--wgs::intervals ${wgs_intervals}" : ''

    """
    riker multi \\
        -i ${bam} \\
        ${reference_arg} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        ${hybcap_opts} \\
        ${error_vcf_arg} \\
        ${error_intervals_arg} \\
        ${gcbias_exclude_intervals_arg} \\
        ${rna_gene_model_arg} \\
        ${rna_ribosomal_intervals_arg} \\
        ${wgs_intervals_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment-metrics.txt
    touch ${prefix}.base-distribution-by-cycle.txt
    touch ${prefix}.mean-quality-by-cycle.txt
    touch ${prefix}.quality-score-distribution.txt
    touch ${prefix}.error-mismatch.txt
    touch ${prefix}.error-overlap.txt
    touch ${prefix}.error-indel.txt
    touch ${prefix}.gcbias-detail.txt
    touch ${prefix}.gcbias-summary.txt
    touch ${prefix}.hybcap-metrics.txt
    touch ${prefix}.hybcap-per-target.txt
    touch ${prefix}.hybcap-per-base.txt
    touch ${prefix}.isize-metrics.txt
    touch ${prefix}.isize-histogram.txt
    touch ${prefix}.wgs-metrics.txt
    touch ${prefix}.wgs-coverage.txt
    touch ${prefix}.base-distribution-by-cycle.pdf
    touch ${prefix}.gcbias-chart.pdf
    touch ${prefix}.isize-histogram.pdf
    touch ${prefix}.mean-quality-by-cycle.pdf
    touch ${prefix}.quality-score-distribution.pdf
    touch ${prefix}.wgs-coverage.pdf
    """
}
