process SENTIEON_HSMETRICS {
    tag "${meta.id}"
    label 'process_low'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7ee64f3b4cd58eaa6ed3a7a0769c0a5f4fbd150842fa64358831970bacaf37e6/data'
        : 'community.wave.seqera.io/library/sentieon:202503.03--df1987151f8b6d33'}"

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
