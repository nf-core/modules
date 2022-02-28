def VERSION = '377' // Version information not provided by tool on CLI

process UCSC_BIGWIGAVERAGEOVERBED {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ucsc-bigwigaverageoverbed=377" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bigwigaverageoverbed:377--h0b8a92a_2' :
        'quay.io/biocontainers/ucsc-bigwigaverageoverbed:377--h0b8a92a_2' }"

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
