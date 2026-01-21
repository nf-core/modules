process FASTK_HISTEX {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b5/b5b07773b60921f43e2839edab692377dcd725d4041adff5520838154e46f487/data' :
        'community.wave.seqera.io/library/fastk:1.2--580652dfcc8e7a12' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("*.hist.txt"), emit: hist
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), eval('echo 1.2'), emit: versions_fastk, topic: versions

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
