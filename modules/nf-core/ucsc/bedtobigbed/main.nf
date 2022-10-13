process UCSC_BEDTOBIGBED {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::ucsc-bedtobigbed=377" : null)
    def container_image = "/ucsc-bedtobigbed:377--ha8a8165_3"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

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
