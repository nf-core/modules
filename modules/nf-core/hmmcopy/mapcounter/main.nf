process HMMCOPY_MAPCOUNTER {
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::hmmcopy=0.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_7':
        'biocontainers/hmmcopy:0.1.1--h2e03b76_7' }"

    input:
    tuple val(meta), path(bigwig)

    output:
    tuple val(meta), path("*.wig"), emit: wig
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mapCounter \\
        $args \\
        $bigwig > ${prefix}.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmcopy: \$(echo $VERSION)
    END_VERSIONS
    """
}
