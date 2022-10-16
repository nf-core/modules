process SEROBA_RUN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::seroba=1.0.2" : null)
    def container_image = "seroba:1.0.2--pyhdfd78af_1"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')

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
