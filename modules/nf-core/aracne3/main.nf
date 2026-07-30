process ARACNE3 {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/papaemmelab/aracne3:v1.0.1"

    input:
    tuple val(meta), path(expression_matrix)
    path regulators

    output:
    tuple val(meta), path("*${meta.id}.tsv"),   emit: consensus_network
    tuple val(meta), path("subnets"),           emit: subnets

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    /ARACNe3/build/src/app/ARACNe3_app_release \\
        -e $expression_matrix \\
        -r $regulators \\
        -o . \\
        --runid $meta.id \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''

    """
    cat <<-END_VERSIONS > versions.yml
    /ARACNe3/build/src/app/ARACNe3_app_release \\
        -e $expression_matrix \\
        -r $regulators \\
        -o . \\
        --runid $meta.id \\
        $args

    END_VERSIONS
    """
}
