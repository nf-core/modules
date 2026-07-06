process RIKER_MULTI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/54/54c7820b49cfb5fada32c1825ac8a46c05d0085105c50055cbce4c37701ab95e/data' :
        'community.wave.seqera.io/library/riker:0.4.0--4e7eeb0beed906c0' }"

    input:
    tuple val(meta),  path(bam), path(bai), path(baits, stageAs: 'baits/*'), path(targets, stageAs: 'targets/*')
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.alignment-metrics.txt"),          emit: alignment_metrics, optional: true
    tuple val(meta), path("*.base-distribution-by-cycle.txt"), emit: base_dist,         optional: true
    tuple val(meta), path("*.mean-quality-by-cycle.txt"),      emit: mean_qual,         optional: true
    tuple val(meta), path("*.quality-score-distribution.txt"), emit: qual_dist,         optional: true
    tuple val(meta), path("*.error-mismatch.txt"),             emit: error_mismatch,    optional: true
    tuple val(meta), path("*.error-overlap.txt"),              emit: error_overlap,     optional: true
    tuple val(meta), path("*.error-indel.txt"),                emit: error_indel,       optional: true
    tuple val(meta), path("*.gcbias-detail.txt"),              emit: gcbias_detail,     optional: true
    tuple val(meta), path("*.gcbias-summary.txt"),             emit: gcbias_summary,    optional: true
    tuple val(meta), path("*.hybcap-metrics.txt"),             emit: hybcap_metrics,    optional: true
    tuple val(meta), path("*.hybcap-per-target.txt"),          emit: hybcap_per_target, optional: true
    tuple val(meta), path("*.hybcap-per-base.txt*"),           emit: hybcap_per_base,   optional: true
    tuple val(meta), path("*.isize-metrics.txt"),              emit: isize_metrics,     optional: true
    tuple val(meta), path("*.isize-histogram.txt"),            emit: isize_histogram,   optional: true
    tuple val(meta), path("*.wgs-metrics.txt"),                emit: wgs_metrics,       optional: true
    tuple val(meta), path("*.wgs-coverage.txt"),               emit: wgs_coverage,      optional: true
    tuple val(meta), path("*.pdf"),                            emit: pdf,               optional: true
    tuple val("${task.process}"), val('riker'), eval("riker --version 2>&1 | sed 's/riker //'") , topic: versions, emit: versions_riker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def ref         = fasta ? "-r ${fasta}" : ''
    if ((baits as Boolean) ^ (targets as Boolean)) {
        error "RIKER_MULTI: both 'baits' and 'targets' must be provided together, or neither"
    }
    def hybcap_opts = (baits && targets) ? "--hybcap::baits ${baits} --hybcap::targets ${targets}" : ''
    """
    riker multi \\
        -i ${bam} \\
        ${ref} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        ${hybcap_opts} \\
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
