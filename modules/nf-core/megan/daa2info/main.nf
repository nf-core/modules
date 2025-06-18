process MEGAN_DAA2INFO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/megan:6.24.20--h9ee0642_0':
        'biocontainers/megan:6.24.20--h9ee0642_0' }"

    input:
    tuple val(meta), path(daa)
    val(megan_summary)

    output:
    tuple val(meta), path("*.txt.gz"), emit: txt_gz
    tuple val(meta), path("*.megan") , emit: megan, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def summary = megan_summary ? "-es ${prefix}.megan" : ""

    """
    daa2info \\
        -i ${daa} \\
        -o ${prefix}.txt.gz \\
        ${summary} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megan: \$(echo \$(rma2info 2>&1) | grep version | sed 's/.*version //g;s/, built.*//g')
    END_VERSIONS
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def megan_cmd = megan_summary ? "touch ${prefix}.megan" : ""

    """
    echo | gzip > ${prefix}.txt.gz
    ${megan_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megan: \$(echo \$(rma2info 2>&1) | grep version | sed 's/.*version //g;s/, built.*//g')
    END_VERSIONS
    """
}
