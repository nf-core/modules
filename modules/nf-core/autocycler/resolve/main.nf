process AUTOCYCLER_RESOLVE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autocycler:0.5.2--h3ab6199_0':
        'biocontainers/autocycler:0.5.2--h3ab6199_0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("resolve/${prefix}/${prefix}_3_bridged.gfa"), emit: bridged
    tuple val(meta), path("resolve/${prefix}/${prefix}_4_merged.gfa"),  emit: merged
    tuple val(meta), path("resolve/${prefix}/${prefix}_5_final.gfa"),   emit: resolved
    tuple val("${task.process}"), val("autocycler"), eval("autocycler --version |  sed 's/^[^ ]* //'"), emit: versions_autocycler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    autocycler resolve \\
        $args \\
        -c .

    mkdir resolve/$prefix -p
    mv 3_bridged.gfa resolve/${prefix}/${prefix}_3_bridged.gfa
    mv 4_merged.gfa  resolve/${prefix}/${prefix}_4_merged.gfa
    mv 5_final.gfa   resolve/${prefix}/${prefix}_5_final.gfa
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p resolve/$prefix
    touch resolve/${prefix}/${prefix}_3_bridged.gfa
    touch resolve/${prefix}/${prefix}_4_merged.gfa
    touch resolve/${prefix}/${prefix}_5_final.gfa
    """
}
