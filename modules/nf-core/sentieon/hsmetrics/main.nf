process SENTIEON_HSMETRICS {
    tag "${meta.id}"
    label 'process_low'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f1dfe59ef66d7326b43db9ab1f39ce6220b358a311078c949a208f9c9815d4e/data'
        : 'community.wave.seqera.io/library/sentieon:202503.01--1863def31ed8e4d5'}"

    input:
    tuple val(meta), path(bam), path(bai), path(bait_intervals), path(target_intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path('*.txt'), emit: metrics
    tuple val("${task.process}"), val('sentieon'), eval("sentieon driver --version | sed 's/.*-//g'"), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = bam.sort().collect {bam_ -> "-i ${bam_}" }.join(' ')
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
        --algo HsMetricAlgo \\
        --targets_list ${target_intervals} \\
        --baits_list ${bait_intervals} \\
        ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
