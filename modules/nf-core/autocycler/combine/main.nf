process AUTOCYCLER_COMBINE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a59a19abe50e65d964ea723d2c310a0005ef5a965beb7bc740d02894dd056e8/data'
:         'community.wave.seqera.io/library/autocycler:0.5.2--91b58d5b8fccdfb9' }"

    input:
    tuple val(meta), path(clusters)

    output:
    tuple val(meta), path("combine/${prefix}/consensus_assembly.fasta"), emit: fasta
    tuple val(meta), path("combine/${prefix}/consensus_assembly.gfa"),   emit: gfa
    tuple val(meta), path("combine/${prefix}/consensus_assembly.yaml"),  emit: stats
    tuple val("${task.process}"), val("autocycler"), eval("autocycler --version |  sed 's/^[^ ]* //'"), emit: versions_autocycler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    autocycler combine \\
        $args \\
        -i $clusters \\
        -a combine

    mkdir combine/$prefix
    mv combine/consensus_assembly.fasta combine/${prefix}/consensus_assembly.fasta
    mv combine/consensus_assembly.gfa combine/${prefix}/consensus_assembly.gfa
    mv combine/consensus_assembly.yaml combine/${prefix}/consensus_assembly.yaml
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p combine/$prefix
    touch combine/${prefix}/consensus_assembly.{fasta,gfa,yaml}
    """
}
