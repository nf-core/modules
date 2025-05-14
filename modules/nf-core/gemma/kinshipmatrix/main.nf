process GEMMA_KINSHIPMATRIX {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gemma:0.98.5--ha36d3ea_0':
        'biocontainers/gemma/0.98.1dfsg-1-deb' }"

    input:
    tuple val(meta), path(genotype)
    tuple val(meta2), path(phenotype)

    output:
    // TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct
    tuple val(meta), path("*.{}"), emit: heat_map
    tuple val(meta), path("*.{}"), emit: quality_control_report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gemma \\
        $args \\
        -@ $task.cpus \\
        -g $genotype \\
        -p $phenotype  \\
        -gk
        -o $meta\.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch $meta\.out
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma --version)
    END_VERSIONS
    """
}
