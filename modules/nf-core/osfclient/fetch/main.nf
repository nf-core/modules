process OSFCLIENT_FETCH {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81145c0eaff8fd8858115f16123f0538a8d5da0a9ad3cbc9672ea82a41f1af5f/data'
        : 'community.wave.seqera.io/library/osfclient:0.0.5--f06e6c22843b8c2a'}"

    input:
    tuple val(meta), val(project_id), val(path)

    output:
    tuple val(meta), path("${outname}"), emit: download_files
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    outname = path.tokenize('/').last()
    """
    osf \\
        -p ${project_id} \\
        fetch ${path} \\
        ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        osfclient: \$(osf --version  | sed 's/osf //g')
    END_VERSIONS
    """

    stub:
    outname = path.tokenize('/').last()
    """
    touch ${outname}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        osfclient: \$(osf --version  | sed 's/osf //g')
    END_VERSIONS
    """
}
