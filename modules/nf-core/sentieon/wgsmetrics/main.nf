process SENTIEON_WGSMETRICS {
    tag "${meta.id}"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(intervals_list)

    output:
    tuple val(meta), path('*.txt'), emit: wgs_metrics
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval = intervals_list ? "--interval ${intervals_list}" : ""
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
        ${interval} \\
        ${args} \\
        --algo WgsMetricsAlgo \\
        ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
