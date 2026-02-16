process NAIL_SEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nail:0.4.0--h4349ce8_0':
        'biocontainers/nail:0.4.0--h4349ce8_0' }"

    input:
    tuple val(meta) , path(query)
    tuple val(meta2), path(target)
    val write_align

    output:
    tuple val(meta), path("${prefix}.txt"), emit: output
    tuple val(meta), path("${prefix}.tbl"), emit: target_summary
    tuple val(meta), path("${prefix}.ali"), emit: alignments, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    alignment = write_align ? "--ali-out ${prefix}.ali" : '' // no def here due to current Nextflow bug; cause: Variable `prefix` already defined in the process scope
    """
    nail search \\
        $args \\
        -t $task.cpus \\
        $alignment \\
        --tbl-out ${prefix}.tbl \\
        ${query} \\
        ${target} > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nail: \$(echo \$(nail --version 2>&1) | sed 's/^.*nail\\w*//' )
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    touch ${prefix}.tbl
    ${write_align ? "touch ${prefix}.ali" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nail: \$(echo \$(nail --version 2>&1) | sed 's/^.*nail\\w*//' )
    END_VERSIONS
    """
}
