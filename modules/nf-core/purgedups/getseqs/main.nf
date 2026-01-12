process PURGEDUPS_GETSEQS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--py39h7132678_1':
        'biocontainers/purge_dups:1.2.6--py39h7132678_1' }"

    input:
    tuple val(meta), path(assembly), path(bed)

    output:
    tuple val(meta), path("*.hap.fa")   , emit: haplotigs
    tuple val(meta), path("*.purged.fa"), emit: purged
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.6' // WARN: Incorrect version printed inside the container, please check this if bumping version
    """
    get_seqs \\
        $args \\
        -e $bed \\
        -p $prefix \\
        $assembly

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.6' // WARN: Incorrect version printed inside the container, please check this if bumping version
    """
    touch "${prefix}.hap.fa"
    touch "${prefix}.purged.fa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: $VERSION
    END_VERSIONS
    """
}
