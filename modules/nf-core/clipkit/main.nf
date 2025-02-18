process CLIPKIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clipkit:2.4.0--pyhdfd78af_0':
        'biocontainers/clipkit:2.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(aln)

    output:
    tuple val(meta), path("*.clipkit"), emit: clipkit
    tuple val(meta), path("log.txt")  , emit: log
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clipkit \\
        $args \\
        $aln \\
        -o ${prefix}.clipkit > log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clipkit: \$(clipkit --version |& sed '1!d ; s/clipkit //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clipkit
    touch log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clipkit: \$(clipkit --version |& sed '1!d ; s/clipkit //')
    END_VERSIONS
    """
}
