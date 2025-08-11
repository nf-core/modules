process PURGEDUPS_CALCUTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--py39h7132678_1':
        'biocontainers/purge_dups:1.2.6--py39h7132678_1' }"

    input:
    tuple val(meta), path(stat)

    output:
    tuple val(meta), path("*.cutoffs")    , emit: cutoff
    tuple val(meta), path("*.calcuts.log"), emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PURGEDUPS modules give segmentation faults when testing using conda and so are currently not recommended"
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.6' // WARN: Incorrect version printed inside the container, please check this if bumping version
    """
    calcuts $args $stat > ${prefix}.cutoffs 2>| >(tee ${prefix}.calcuts.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: $VERSION
    END_VERSIONS
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PURGEDUPS modules give segmentation faults when testing using conda and so are currently not recommended"
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.6' // WARN: Incorrect version printed inside the container, please check this if bumping version ( \$( purge_dups -h |& sed '3!d; s/.*: //' ))
    """
    touch "${prefix}.cutoffs"
    touch "${prefix}.calcuts.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: $VERSION
    END_VERSIONS
    """
}
