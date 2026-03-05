process SENTIEON_DATAMETRICS {
    tag "${meta.id}"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/ae67a134620c3af22c8563a2913c4639caa0d75ce25764e7b10c996b242aa023/data'
        : 'community.wave.seqera.io/library/sentieon_gnuplot:41931fca35668c97'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    val plot_results

    output:
    tuple val(meta), path('*mq_metrics.txt'),  emit: mq_metrics
    tuple val(meta), path('*qd_metrics.txt'),  emit: qd_metrics
    tuple val(meta), path('*gc_summary.txt'),  emit: gc_summary
    tuple val(meta), path('*gc_metrics.txt'),  emit: gc_metrics
    tuple val(meta), path('*aln_metrics.txt'), emit: aln_metrics
    tuple val(meta), path('*is_metrics.txt'),  emit: is_metrics
    tuple val(meta), path('*mq_metrics.pdf'),  emit: mq_plot, optional: true
    tuple val(meta), path('*qd_metrics.pdf'),  emit: qd_plot, optional: true
    tuple val(meta), path('*is_metrics.pdf'),  emit: is_plot, optional: true
    tuple val(meta), path('*gc_metrics.pdf'),  emit: gc_plot, optional: true
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = bam.sort().collect {in -> "-i ${in}" }.join(' ')
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon \\
        driver \\
        -t ${task.cpus} \\
        -r ${fasta} \\
        ${input} \\
        ${args} \\
        --algo GCBias --summary ${prefix}_gc_summary.txt ${prefix}_gc_metrics.txt \\
        --algo MeanQualityByCycle ${prefix}_mq_metrics.txt \\
        --algo QualDistribution ${prefix}_qd_metrics.txt \\
        --algo InsertSizeMetricAlgo ${prefix}_is_metrics.txt  \\
        --algo AlignmentStat ${prefix}_aln_metrics.txt

    if ${plot_results}
    then
        sentieon plot GCBias -o ${prefix}_gc_metrics.pdf ${prefix}_gc_metrics.txt
        sentieon plot MeanQualityByCycle -o ${prefix}_mq_metrics.pdf ${prefix}_mq_metrics.txt
        sentieon plot QualDistribution -o ${prefix}_qd_metrics.pdf  ${prefix}_qd_metrics.txt
        sentieon plot InsertSizeMetricAlgo -o ${prefix}_is_metrics.pdf ${prefix}_is_metrics.txt
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mq_metrics.txt
    touch ${prefix}_qd_metrics.txt
    touch ${prefix}_gc_summary.txt
    touch ${prefix}_gc_metrics.txt
    touch ${prefix}_aln_metrics.txt
    touch ${prefix}_is_metrics.txt

    if ${plot_results}
    then
        touch ${prefix}_gc_metrics.pdf
        touch ${prefix}_mq_metrics.pdf
        touch ${prefix}_qd_metrics.pdf
        touch ${prefix}_is_metrics.pdf
    fi
    """
}
