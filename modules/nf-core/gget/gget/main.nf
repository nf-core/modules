process GGET_GGET {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gget=0.27.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gget:0.27.2--pyhdfd78af_0':
        'quay.io/biocontainers/gget:0.27.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*")     , emit: file, optional: true
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val(meta), path("*.csv"), emit: csv,  optional: true
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def inputs = files ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "json"
    """
    gget \\
        $args \\
        -o ${prefix}.${suffix} \\
        $inputs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gget: \$(echo \$(gget --version 2>&1 | sed 's/gget version: //g'))
    END_VERSIONS
    """
}
