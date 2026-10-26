process RATTLE_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/597c4110442b53dc0d73420af28bc774518719d974b2be3190dfa478da17acac/data'
        : 'community.wave.seqera.io/library/rattle:1.0--ae3dd63e7a98745b'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("clusters.out"), emit: clusters
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('rattle'), val("1.0"), emit: versions_rattle, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rattle \\
        cluster ${args} \\
        -t ${task.cpus} \\
        -i ${reads}
    """

    stub:
    """
    touch clusters.out
    """
}
