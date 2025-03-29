process GGET_GGET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gget:0.29.0--pyhdfd78af_0':
        'biocontainers/gget:0.29.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*[!versions.yml][!${prefix}.${extension}]*"), emit: files, optional: true
    tuple val(meta), path("${prefix}.${extension}")                   , emit: output, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def inputs = files ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.startsWith("pdb")  ? "pdb" :
                args.contains("-csv")   ? "csv" :
                "json"
    """
    gget \\
        $args \\
        -o ${prefix}.${extension} \\
        $inputs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gget: \$(echo \$(gget --version 2>&1 | sed 's/gget version: //g'))
    END_VERSIONS
    """
}
