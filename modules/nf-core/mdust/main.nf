process MDUST {
    tag "$meta.id"
    label 'process_single'

    // WARN: Manually update when changing Bioconda assets
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mdust:2006.10.17--h470a237_1':
        'biocontainers/mdust:2006.10.17--h470a237_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_mdust"
    def VERSION = '2006.10.17' // WARN: Manually update when changing Bioconda assets
    if( "$fasta" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    mdust \\
        $fasta \\
        $args \\
        > "${prefix}.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mdust: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_mdust"
    def VERSION = '2006.10.17' // WARN: Manually update when changing Bioconda assets
    if( "$fasta" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mdust: $VERSION
    END_VERSIONS
    """
}
