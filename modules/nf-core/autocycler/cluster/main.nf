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
    tuple val(meta), path("clustering/$prefix/qc_pass/*/*.gfa"),  emit: clusters
    tuple val(meta), path("clustering/$prefix/qc_pass/*/*.yaml"), emit: clusterstats
    tuple val(meta), path("clustering/$prefix/*.newick"),         emit: newick
    tuple val(meta), path("clustering/$prefix/*.tsv"),            emit: tsv
    tuple val(meta), path("clustering/$prefix/*.phylip"),         emit: pairwisedistances
    tuple val(meta), path("clustering/$prefix/*.yaml"),           emit: stats
    tuple val("${task.process}"), val("autocycler"), eval("autocycler --version |  sed 's/^[^ ]* //'"), emit: versions_autocycler, topic: versions

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
    mv clustering/* $prefix
    mv $prefix clustering
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    mkdir clustering/$prefix/qc_pass/cluster_000 -p
    touch clustering/$prefix/clustering.{newick,yaml,tsv}
    touch clustering/$prefix/pairwise_distances.phylip
    touch clustering/$prefix/qc_pass/cluster_000/0_untrimmed.{gfa,yaml}
    """
}
