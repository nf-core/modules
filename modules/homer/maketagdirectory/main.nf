def VERSION = '4.11' // Version information not provided by tool on CLI

process HOMER_MAKETAGDIRECTORY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::homer=4.11=pl526hc9558a2_3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(bed)
    path fasta

    output:
    tuple val(meta), path("tag_dir"), emit: tagdir
    path  "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    makeTagDirectory \\
        tag_dir \\
        $args \\
        $bed \\
        -genome $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
