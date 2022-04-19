def VERSION = '2.1' // Version information not provided by tool on CLI

process GAMMA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gamma=2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gamma%3A2.1--hdfd78af_0':
        'quay.io/biocontainers/gamma:2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("*.gamma"), emit: gamma
    tuple val(meta), path("*.psl")  , emit: psl
    tuple val(meta), path("*.gff")  , emit: gff      , optional:true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    GAMMA.py \\
    $args \\
    $fasta \\
    $db \\
    $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gamma: $VERSION
    END_VERSIONS
    """
}
