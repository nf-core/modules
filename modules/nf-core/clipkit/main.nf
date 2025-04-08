process CLIPKIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clipkit:2.4.1--pyhdfd78af_0':
        'biocontainers/clipkit:2.4.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(aln)
    val out_format

    output:
    tuple val(meta), path("${prefix}.${out_extension}"), emit: clipkit
    tuple val(meta), path("${prefix}.log")             , emit: log
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    out_extension = out_format ? out_format : "clipkit"
    """
    clipkit \\
        $args \\
        $aln \\
        -o ${prefix}.${out_extension} > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clipkit: \$(clipkit --version |& sed '1!d ; s/clipkit //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    out_extension = out_format ? out_format : "clipkit"
    """
    touch ${prefix}.${out_extension}
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clipkit: \$(clipkit --version |& sed '1!d ; s/clipkit //')
    END_VERSIONS
    """
}
