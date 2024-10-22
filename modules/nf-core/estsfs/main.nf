process ESTSFS {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/est-sfs:2.04--h245ed52_0':
        'biocontainers/est-sfs:2.04--h245ed52_0' }"

    input:
    tuple val(meta), path(e_config), path(data), path(seed)

    output:
    tuple val(meta), path("${prefix}_sfs.txt")   , emit: sfs_out
    tuple val(meta), path("${prefix}_pvalues.txt"), emit: pvalues_out
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '2.04' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    est-sfs ${e_config} ${data} ${seed} ${prefix}_sfs.txt ${prefix}_pvalues.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        est-sfs: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.04' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_sfs.txt
    touch ${prefix}_pvalues.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        est-sfs: $VERSION
    END_VERSIONS
    """
}
