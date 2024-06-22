process MD5SUM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*.md5"), emit: checksum
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    IFS=\$(echo -en "\n\b")
    for FILE in $files
    do

    md5sum \\
        $args \\
        \${FILE} \\
        > "\${FILE}.md5"

    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        md5sum: \$(echo \$(md5sum --version 2>&1 | head -n 1| sed 's/^.*) //;' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    IFS=\$(echo -en "\n\b")
    for FILE in $files
    do
    touch "\${FILE}.md5"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        md5sum: \$(echo \$(md5sum --version 2>&1 | head -n 1| sed 's/^.*) //;' ))
    END_VERSIONS
    """

}
