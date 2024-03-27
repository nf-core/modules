process XENGSORT_INDEX {
    tag "$host_fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xengsort:2.0.5--pyhdfd78af_0':
        'biocontainers/xengsort:2.0.5--pyhdfd78af_0' }"

    input:
    path(host_fasta, stageAs: "host/*")
    path(graft_fasta, stageAs: "graft/*")
    val index
    val nobjects
    val mask

    output:
    path "${index}.hash"          , emit: hash
    path "${index}.info"          , emit: info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    xengsort \\
        index \\
        $args \\
        --index $index \\
        --host $host_fasta \\
        --graft $graft_fasta \\
        --nobjects $nobjects \\
        --mask '$mask' \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xengsort: \$(xengsort --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${index}.info
    touch ${index}.hash

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xengsort: \$(xengsort --version)
    END_VERSIONS
    """
}
