process WGSIM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::wgsim=1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wgsim:1.0--h5bf99c6_4':
        'biocontainers/wgsim:1.0--h5bf99c6_4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fastq"), emit: fastq
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def WGSIM_VERSION = "0.3.1-r13" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    wgsim \\
        $args \\
        $fasta \\
        ${prefix}_1.fastq ${prefix}_2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wgsim: $WGSIM_VERSION
    END_VERSIONS
    """
}
