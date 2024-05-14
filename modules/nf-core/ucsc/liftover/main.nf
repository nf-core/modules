process UCSC_LIFTOVER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::ucsc-liftover=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-liftover:377--h0b8a92a_3' :
        'biocontainers/ucsc-liftover:377--h0b8a92a_3' }"

    input:
    tuple val(meta), path(bed)
    path(chain)

    output:
    tuple val(meta), path("*.lifted.bed")  , emit: lifted
    tuple val(meta), path("*.unlifted.bed"), emit: unlifted
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    liftOver \\
        $args \
        $bed \\
        $chain \\
        ${prefix}.lifted.bed \\
        ${prefix}.unlifted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
