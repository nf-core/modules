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
    tuple val(meta), path("$prefix/clustering/qc_pass/*/*.gfa"),  emit: clusters
    tuple val(meta), path("$prefix/clustering/qc_pass/*/*.yaml"), emit: clusterstats
    tuple val(meta), path("$prefix/clustering/*.newick"),         emit: newick
    tuple val(meta), path("$prefix/clustering/*.tsv"),            emit: tsv
    tuple val(meta), path("$prefix/clustering/*.phylip"),         emit: pairwisedistances
    tuple val(meta), path("$prefix/clustering/*.yaml"),           emit: stats
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    autocycler cluster \\
        $args \\
        -a .

    mkdir $prefix
    mv clustering $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    mkdir $prefix/clustering/qc_pass/cluster_000 -p
    touch $prefix/clustering/clustering.{newick,yaml,tsv}
    touch $prefix/clustering/pairwise_distances.phylip
    touch $prefix/clustering/qc_pass/cluster_000/0_untrimmed.{gfa,yaml}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """
}
