process HLALA_PREPAREGRAPH {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::hla-la=1.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hla-la:1.0.3--hd03093a_0':
        'quay.io/biocontainers/hla-la:1.0.3--hd03093a_0' }"

    input:
        tuple val(meta), path(graph)


    output:
    tuple val(meta), path("${graph}")        , emit: graph
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //def zipped = zipped_graph.toString().endsWith(".zip")
    //def unzipped_graph = zipped_graph ? zipped_graph.toString() - ~/\.zip$/: ""


    // OBS: the "../bin/HLA-LA" described in the documentation is found in the docker container in "/usr/local/opt/hla-la/bin/HLA-LA"
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
