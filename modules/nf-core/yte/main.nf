process YTE {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/yte:1.5.4--2466971526dcd2a8':
        'community.wave.seqera.io/library/yte:1.5.4--36d16ca4bab836c1' }"

    input:
    tuple val(meta), path(template)

    output:
    tuple val(meta), path("*.yaml"), emit: rendered
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "1.5.4" // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    yte < ${template} > ${prefix}.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yte: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "1.5.4" // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    touch ${prefix}.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yte: $VERSION
    END_VERSIONS
    """
}
