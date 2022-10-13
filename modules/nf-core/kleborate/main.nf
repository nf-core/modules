process KLEBORATE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kleborate=2.1.0" : null)
    def container_image = "/kleborate:2.1.0--pyhdfd78af_1"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kleborate \\
        $args \\
        --outfile ${prefix}.results.txt \\
        --assemblies $fastas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """
}
