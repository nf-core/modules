process UCSC_BIGWIGAVERAGEOVERBED {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::ucsc-bigwigaverageoverbed=377" : null)
    def container_image = "/ucsc-bigwigaverageoverbed:377--h0b8a92a_2"
                                                            container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bed)
    path bigwig

    output:
    tuple val(meta), path("*.tab"), emit: tab
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    // BUG: bigWigAverageOverBed cannot handle ensembl seqlevels style
    """
    bigWigAverageOverBed \\
        $args \\
        $bigwig \\
        $bed \\
        ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
