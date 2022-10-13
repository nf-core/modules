process KALLISTO_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
        def container_image = "/kallisto:0.46.2--h4f7b962_1"
                                                 container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    path fasta

    output:
    path "kallisto" , emit: idx
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    kallisto \\
        index \\
        $args \\
        -i kallisto \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch kallisto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
    END_VERSIONS
    """
}
