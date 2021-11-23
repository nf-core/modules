process ROARY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::roary=3.13.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/roary:3.13.0--pl526h516909a_0' :
        'quay.io/biocontainers/roary:3.13.0--pl526h516909a_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                    , emit: results
    tuple val(meta), path("results/*.aln"), optional: true, emit: aln
    path "versions.yml"                                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    roary \\
        $args \\
        -p $task.cpus \\
        -f results/ \\
        $gff

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        roary: \$( roary --version )
    END_VERSIONS
    """
}
