process DEEPTOOLS_PLOTPROFILE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0"
    } else {
        container "quay.io/biocontainers/deeptools:3.5.1--py_0"
    }

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tab"), emit: table
    path  "versions.yml"          , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    plotProfile \\
        $options.args \\
        --matrixFile $matrix \\
        --outFileName ${prefix}.plotProfile.pdf \\
        --outFileNameData ${prefix}.plotProfile.tab

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(plotProfile --version | sed -e "s/plotProfile //g")
    END_VERSIONS
    """
}
