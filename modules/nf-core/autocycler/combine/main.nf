process AUTOCYCLER_COMBINE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autocycler:0.5.2--h3ab6199_0':
        'biocontainers/autocycler:0.5.2--h3ab6199_0' }"

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
