process AUTOCYCLER_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autocycler:0.5.2--h3ab6199_0':
        'biocontainers/autocycler:0.5.2--h3ab6199_0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("clustering/qc_pass/*/*.gfa"),  emit: clusters
    tuple val(meta), path("clustering/qc_pass/*/*.yaml"), emit: clusterstats
    tuple val(meta), path("clustering/*.newick"),         emit: newick
    tuple val(meta), path("clustering/*.tsv"),            emit: tsv
    tuple val(meta), path("clustering/*.phylip"),         emit: pairwisedistances
    tuple val(meta), path("clustering/*.yaml"),           emit: stats
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    """
    autocycler cluster \\
        $args \\
        -a .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    """
    mkdir clustering/qc_pass/cluster_000 -p
    touch clustering/clustering.{newick,yaml,tsv}
    touch clustering/pairwise_distances.phylip
    touch clustering/qc_pass/cluster_000/0_untrimmed.{gfa,yaml}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """
}
