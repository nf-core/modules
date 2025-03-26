process PAML_BUILDCTL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data':
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"


    input:
    tuple val(meta), path(phy)
    tuple val(meta), path(tree)

    output:
    tuple val(meta), path("settings.ctl")   , emit: ctl
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "seqfile = ${phy}" > settings.ctl
    echo "treefile = ${tree}" >> settings.ctl
    echo "outfile = result.txt" >> settings.ctl
    echo "$args" >> settings.ctl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch settings.ctl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version)
    END_VERSIONS
    """
}
