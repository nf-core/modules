process GAPPA_EXAMINEGRAFT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gappa=0.8.0" : null)
    def container_image = "gappa:0.8.0--h9a82719_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(jplace)

    output:
    tuple val(meta), path("*.newick"), emit: newick
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gappa \\
        examine \\
        graft \\
        $args \\
        --threads $task.cpus \\
        --file-prefix ${prefix}. \\
        --jplace-path $jplace

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gappa: \$(echo \$(gappa --version 2>&1 | sed 's/v//' ))
    END_VERSIONS
    """
}
