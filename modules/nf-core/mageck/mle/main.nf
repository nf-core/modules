process MAGECK_MLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageck:0.5.9.5--py39h1f90b4d_3':
        'biocontainers/mageck:0.5.9.5--py39h1f90b4d_3' }"

    input:
    tuple val(meta), path(count_table)
    path(design_matrix)

    output:
    tuple val(meta), path("*.gene_summary.txt") , emit: gene_summary
    tuple val(meta), path("*.sgrna_summary.txt"), emit: sgrna_summary
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mageck \\
        mle \\
        $args \\
        --threads $task.cpus \\
        -k $count_table \\
        -d $design_matrix \\
        -n $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.gene_summary.txt
    touch ${prefix}.sgrna_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """


}
