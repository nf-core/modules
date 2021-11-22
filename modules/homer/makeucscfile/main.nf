def VERSION = '4.11'

process HOMER_MAKEUCSCFILE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::homer=4.11=pl526hc9558a2_3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(tagDir)

    output:
    tuple val(meta), path("tag_dir/*ucsc.bedGraph.gz"), emit: bedGraph
    path  "versions.yml"                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    makeUCSCfile \\
        $tagDir \\
        -o auto
        $args

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        homer: \$(echo $VERSION)
    END_VERSIONS
    """
}
