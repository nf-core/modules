process PAIRTOOLS_RESTRICT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::pairtools=0.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pairtools:0.3.0--py37hb9c2fc3_5"
    } else {
        container "quay.io/biocontainers/pairtools:0.3.0--py37hb9c2fc3_5"
    }

    input:
    tuple val(meta), path(pairs)
    path frag

    output:
    tuple val(meta), path("*.pairs.gz"), emit: restrict
    path "versions.yml"                , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    pairtools \\
        restrict \\
        -f $frag \\
        $options.args \\
        -o ${prefix}.pairs.gz \\
        $pairs

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
    END_VERSIONS
    """
}
