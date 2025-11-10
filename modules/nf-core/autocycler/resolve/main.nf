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
    tuple val(meta), path("$prefix/3_bridged.gfa"), emit: bridged
    tuple val(meta), path("$prefix/4_merged.gfa"),  emit: merged
    tuple val(meta), path("$prefix/5_final.gfa"),   emit: resolved
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    autocycler resolve \\
        $args \\
        -c .

    mkdir $prefix
    mv 3_bridged.gfa $prefix
    mv 4_merged.gfa $prefix
    mv 5_final.gfa $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix
    touch $prefix/3_bridged.gfa
    touch $prefix/4_merged.gfa
    touch $prefix/5_final.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """
}
