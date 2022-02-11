def VERSION = '377' // Version information not provided by tool on CLI

process UCSC_WIGTOBIGWIG {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ucsc-wigtobigwig=377" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-wigtobigwig:377--h0b8a92a_2' :
        'quay.io/biocontainers/ucsc-wigtobigwig:377--h0b8a92a_2' }"

    input:
    tuple val(meta), path(wig)
    path sizes

    output:
    tuple val(meta), path("*.bw"), emit: bw
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    wigToBigWig \\
        $args \\
        $wig \\
        $sizes \\
        ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
