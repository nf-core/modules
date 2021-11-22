process PAIRTOOLS_FLIP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pairtools=0.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:0.3.0--py37hb9c2fc3_5' :
        'quay.io/biocontainers/pairtools:0.3.0--py37hb9c2fc3_5' }"

    input:
    tuple val(meta), path(sam)
    path chromsizes

    output:
    tuple val(meta), path("*.flip.gz"), emit: flip
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    pairtools \\
        flip \\
        -c $chromsizes \\
        $args \\
        -o ${prefix}.flip.gz \\
        $sam

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        ${getSoftwareName(task.process)}: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
    END_VERSIONS
    """
}
