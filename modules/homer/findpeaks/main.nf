def VERSION = '4.11' // Version information not provided by tool on CLI

process HOMER_FINDPEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::homer=4.11=pl526hc9558a2_3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(tagDir)

    output:
    tuple val(meta), path("*peaks.txt"), emit: txt
    path  "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """

    findPeaks \\
        $tagDir \\
        $args \\
        -o ${prefix}.peaks.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        homer: $VERSION
    END_VERSIONS
    """
}
