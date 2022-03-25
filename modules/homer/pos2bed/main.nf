def VERSION = '4.11' // Version information not provided by tool on CLI

process HOMER_POS2BED {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::homer=4.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(peaks)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pos2bed.pl $peaks > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
