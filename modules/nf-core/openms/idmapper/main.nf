process OPENMS_IDMAPPER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.1.0--h8964181_3':
        'biocontainers/openms:3.1.0--h8964181_3' }"

    input:
    tuple val(meta), path(id_file), path(map_file)

    output:
    path "*.${id_file.getExtension()}", emit: xml
    path "versions.yml", emit: version
    path "*.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    IDMapper \\
        -id ${id_file} \\
        -in ${map_file} \\
        -threads $task.cpus \\
        -out ${prefix}.${id_file.getExtension()} \\
        $args \\
        2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IDMapper: \$(IDMapper 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${id_file.getExtension()}
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IDMapper: \$(IDMapper 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
    END_VERSIONS
    """
}
