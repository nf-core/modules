process UCSC_BEDTOBIGBED {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::ucsc-bedtobigbed=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:377--ha8a8165_3' :
        'biocontainers/ucsc-bedtobigbed:377--ha8a8165_3' }"

    input:
    tuple val(meta), path(bed)
    path  sizes
    path  autosql

    output:
    tuple val(meta), path("*.bigBed"), emit: bigbed
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def as_option = autosql ? "-as=${autosql}" : ""
    def VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bedToBigBed \\
        $bed \\
        $sizes \\
        $as_option \\
        $args \\
        ${prefix}.bigBed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
