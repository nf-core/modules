process SEROBA_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seroba=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seroba:1.0.2--pyhdfd78af_1':
        'biocontainers/seroba:1.0.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}/${prefix}.tsv")                              , emit: tsv
    tuple val(meta), path("${prefix}/detailed_serogroup_info.txt"), optional: true, emit: txt
    path "versions.yml"                                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    seroba \\
        runSerotyping \\
        $reads \\
        $prefix \\
        $args

    # Avoid name collisions
    mv ${prefix}/pred.tsv ${prefix}/${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seroba: \$(seroba version)
    END_VERSIONS
    """
}
