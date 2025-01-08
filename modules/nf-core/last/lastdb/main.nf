process LAST_LASTDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/last:1595--h43eeafb_0' :
        'biocontainers/last:1595--h43eeafb_0' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("lastdb"), emit: index
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir lastdb
    lastdb \\
        $args \\
        -P $task.cpus \\
        lastdb/${prefix} \\
        $fastx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version 2>&1 | sed 's/lastdb //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir lastdb
    touch lastdb/${prefix}.bck
    touch lastdb/${prefix}.des
    touch lastdb/${prefix}.prj
    touch lastdb/${prefix}.sds
    touch lastdb/${prefix}.ssp
    touch lastdb/${prefix}.suf
    touch lastdb/${prefix}.tis

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version 2>&1 | sed 's/lastdb //')
    END_VERSIONS
    """
}
