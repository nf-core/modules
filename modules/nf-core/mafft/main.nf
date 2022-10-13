process MAFFT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::mafft=7.490" : null)
        'https://depot.galaxyproject.org/singularity/mafft:7.490--h779adbc_0':
        "${params.docker_registry ?: 'quay.io/biocontainers'}/mafft:7.490--h779adbc_0" }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fas"), emit: fas
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mafft \\
        --thread ${task.cpus} \\
        ${args} \\
        ${fasta} \\
        > ${prefix}.fas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
