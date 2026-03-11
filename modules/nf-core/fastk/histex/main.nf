process FASTK_HISTEX {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02c05b2ec421debc83883ef9a211291e3220546c12f7c54cb78e66209cb2797d/data' :
        'community.wave.seqera.io/library/fastk:1.2--4bc70c6cd0d420bd' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("*.hist.txt"), emit: hist
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), val('1.2'), emit: versions_fastk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    Histex \\
        $args \\
        $histogram \\
        > ${prefix}.hist.txt
    """

    stub:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hist.txt
    """
}
