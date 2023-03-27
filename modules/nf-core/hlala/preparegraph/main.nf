process HLALA_PREPAREGRAPH {
    tag "$meta.id"
    label 'process_high'



    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hla-la:1.0.3--hd03093a_0':
        'quay.io/biocontainers/hla-la:1.0.3--hd03093a_0' }"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "HLALA_PREPAREGRAPH module does not support Conda. Please use Docker or Singularity."
    }

    input:
    tuple val(meta), path(graph)


    output:
    tuple val(meta), path("${graph}")        , emit: graph
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    /usr/local/opt/hla-la/bin/HLA-LA \\
        --action prepareGraph \\
        --PRG_graph_dir $graph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlala: 1.0.3
    END_VERSIONS
    """
}
