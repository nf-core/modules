process SHIGATYPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shigatyper:2.0.5--pyhdfd78af_0':
        'biocontainers/shigatyper:2.0.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.tsv")     , emit: tsv
    tuple val(meta), path("${prefix}-hits.tsv"), optional: true, emit: hits
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.is_ont) {
        """
        shigatyper \\
            $args \\
            --SE $reads \\
            --ont \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(shigatyper --version | sed -n 's/ShigaTyper \\(.*\\)/\\1/p')
        END_VERSIONS
        """
    } else if (meta.single_end) {
        """
        shigatyper \\
            $args \\
            --SE $reads \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(shigatyper --version | sed -n 's/ShigaTyper \\(.*\\)/\\1/p')
        END_VERSIONS
        """
    } else {
        """
        shigatyper \\
            $args \\
            --R1 ${reads[0]} \\
            --R2 ${reads[1]} \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(shigatyper --version | sed -n 's/ShigaTyper \\(.*\\)/\\1/p')
        END_VERSIONS
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigatyper: \$(shigatyper --version | sed -n 's/ShigaTyper \\(.*\\)/\\1/p')
    END_VERSIONS
    """
}
