// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process AUTOCYCLER_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a59a19abe50e65d964ea723d2c310a0005ef5a965beb7bc740d02894dd056e8/data'
:         'community.wave.seqera.io/library/autocycler:0.5.2--91b58d5b8fccdfb9' }"

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
