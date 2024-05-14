process SHIGATYPER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::shigatyper=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shigatyper%3A2.0.1--pyhdfd78af_0':
        'biocontainers/shigatyper:2.0.1--pyhdfd78af_0' }"

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
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
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
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
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
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    }
}
