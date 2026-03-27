process CRABZ_COMPRESS {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/798b49da13cf3eddbf611aecbcfef584183d80aaec016406a10990ac741d6f8f/data':
        'community.wave.seqera.io/library/crabz:0.10.0--13b570daf913bcc8' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.gz"), emit: archive
    tuple val("${task.process}"), val('crabz'), val('0.10.0'), emit: versions_crabz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${file}"
    """
    crabz \\
        ${args} \\
        -p ${task.cpus} \\
        -o ${prefix}.gz \\
        ${file}
    """

    stub:
    def prefix = task.ext.prefix ?: "${file}"
    """
    echo "" | gzip > ${prefix}.gz
    """
}
