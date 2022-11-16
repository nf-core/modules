process HLALA_PREPAREGRAPH {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::hla-la=1.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hla-la:1.0.3--hd03093a_0':
        'quay.io/biocontainers/hla-la:1.0.3--hd03093a_0' }"

    input:
    path(graph)

    output:
    path(graph)        , emit: folder
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    ../bin/HLA-LA \\
        --action prepareGraph \\
        --PRG_graph_dir $graph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlala: 1.0.3
    END_VERSIONS
    """
}
