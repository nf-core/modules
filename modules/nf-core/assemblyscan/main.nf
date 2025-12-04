process ASSEMBLYSCAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/assembly-scan:1.0.0--pyhdfd78af_0' :
        'biocontainers/assembly-scan:1.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv"),     emit: tsv, optional: true
    tuple val(meta), path("*.json"),    emit: json, optional: true
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def file_format = (task.ext.args?.contains("--json")) ? "json" : "tsv"
    """
    assembly-scan \\
         ${args} \\
         ${assembly} > ${prefix}.${file_format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assemblyscan: \$( assembly-scan --version 2>&1 | sed 's/^.*assembly-scan //; s/Using.*\$//' )
    END_VERSIONS
    """
}
