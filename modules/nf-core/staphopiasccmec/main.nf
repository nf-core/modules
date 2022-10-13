process STAPHOPIASCCMEC {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::staphopia-sccmec=1.0.0" : null)
    def container_image = "/staphopia-sccmec:1.0.0--hdfd78af_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }
n

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    staphopia-sccmec --assembly $fasta $args > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staphopiasccmec: \$(staphopia-sccmec --version 2>&1 | sed 's/^.*staphopia-sccmec //')
    END_VERSIONS
    """
}
