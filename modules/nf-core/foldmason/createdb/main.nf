process FOLDMASON_CREATEDB {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'oras://community.wave.seqera.io/library/foldmason:2.7bd21ed--a45f76ed12b391e6':
            'community.wave.seqera.io/library/foldmason:2.7bd21ed--e7f739473ad6578d' }"

    input:
    tuple val(meta) , path(structures)

    output:
    tuple val(meta), path("${prefix}*"), emit: db
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    foldmason createdb \\
        ${structures} \\
        ${prefix} \\
        $args \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldmason: \$(foldmason | grep "foldmason Version:" | cut -d":" -f 2 | awk '{\$1=\$1;print}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldmason: \$(foldmason | grep "foldmason Version:" | cut -d":" -f 2 | awk '{\$1=\$1;print}')
    END_VERSIONS
    """
}
