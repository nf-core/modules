process GEMMA_KINSHIPMATRIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gemma:0.98.5--ha36d3ea_0':
        'biocontainers/gemma/0.98.5--ha36d3ea_0' }"

    input:
    tuple val(meta), path(genotype)
    tuple val(meta2), path(phenotype)

    output:
    tuple val(meta), path("output/${meta.id}.out.cXX.txt"), emit: matrix
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gemma \\
        $args \\
        -g $genotype \\
        -p $phenotype  \\
        -gk \\
        -o ${meta.id}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir output
    touch output/${meta.id}.out.cXX.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma --version)
    END_VERSIONS
    """
}
