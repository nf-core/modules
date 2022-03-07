process PAIRTOOLS_RESTRICT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::pairtools=0.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:0.3.0--py37hb9c2fc3_5' :
        'quay.io/biocontainers/pairtools:0.3.0--py37hb9c2fc3_5' }"

    input:
    tuple val(meta), path(pairs)
    path frag

    output:
    tuple val(meta), path("*.pairs.gz"), emit: restrict
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pairtools \\
        restrict \\
        -f $frag \\
        $args \\
        -o ${prefix}.pairs.gz \\
        $pairs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
    END_VERSIONS
    """
}
