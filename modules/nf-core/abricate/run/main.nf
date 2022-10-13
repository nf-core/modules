process ABRICATE_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::abricate=1.0.1" : null)
    def container_image = "/abricate%3A1.0.1--ha8f3691_1"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.txt"), emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    abricate \\
        $assembly \\
        $args \\
        --threads $task.cpus > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
