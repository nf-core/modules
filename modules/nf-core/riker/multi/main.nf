process RIKER_MULTI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5bf9ec40db8ba058b6ff37a94ea398f37b766858b3e584016e93643f7dde9f63/data' :
        'community.wave.seqera.io/library/riker:0.2.0--20857cea9478b433' }"

    input:
    tuple val(meta),  path(bam), path(bai), path(baits), path(targets)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.alignment-metrics.txt")             , optional: true, emit: alignment_metrics
    tuple val(meta), path("*.base-distribution-by-cycle.txt")    , optional: true, emit: base_dist
    tuple val(meta), path("*.mean-quality-by-cycle.txt")         , optional: true, emit: mean_qual
    tuple val(meta), path("*.quality-score-distribution.txt")    , optional: true, emit: qual_dist
    tuple val(meta), path("*.error-mismatch.txt")                , optional: true, emit: error_mismatch
    tuple val(meta), path("*.error-overlap.txt")                 , optional: true, emit: error_overlap
    tuple val(meta), path("*.error-indel.txt")                   , optional: true, emit: error_indel
    tuple val(meta), path("*.gcbias-detail.txt")                 , optional: true, emit: gcbias_detail
    tuple val(meta), path("*.gcbias-summary.txt")                , optional: true, emit: gcbias_summary
    tuple val(meta), path("*.hybcap-metrics.txt")                , optional: true, emit: hybcap_metrics
    tuple val(meta), path("*.hybcap-per-target.txt")             , optional: true, emit: hybcap_per_target
    tuple val(meta), path("*.hybcap-per-base.txt*")              , optional: true, emit: hybcap_per_base
    tuple val(meta), path("*.isize-metrics.txt")                 , optional: true, emit: isize_metrics
    tuple val(meta), path("*.isize-histogram.txt")               , optional: true, emit: isize_histogram
    tuple val(meta), path("*.wgs-metrics.txt")                   , optional: true, emit: wgs_metrics
    tuple val(meta), path("*.wgs-coverage.txt")                  , optional: true, emit: wgs_coverage
    tuple val(meta), path("*.pdf")                               , optional: true, emit: pdf
    tuple val("${task.process}"), val('riker'), eval("riker --version 2>&1 | sed 's/riker //'") , topic: versions, emit: versions_riker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def ref         = fasta ? "-r ${fasta}" : ''
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
