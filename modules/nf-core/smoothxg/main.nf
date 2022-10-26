def VERSION = '0.6.1' // Version information not provided by tool on CLI

process SMOOTHXG {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::smoothxg=0.6.1' : null)

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smoothxg:0.6.1--h735f0e5_0' :
        'quay.io/biocontainers/smoothxg:0.6.1--h735f0e5_0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    smoothxg \\
        --threads=$task.cpus \\
        --gfa-in=${gfa} \\
        --smoothed-out=${prefix}.smoothxg.gfa
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smoothxg: $VERSION
    END_VERSIONS
    """
}
