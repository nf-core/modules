process FASTANI {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fastani=1.32" : null)
        def container_image = "/fastani:1.32--he1c1bb9_0"
                                          container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(query)
    path reference

    output:
    tuple val(meta), path("*.ani.txt"), emit: ani
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.batch_input) {
        """
        fastANI \\
            -ql $query \\
            -rl $reference \\
            -o ${prefix}.ani.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
        END_VERSIONS
        """
    } else {
        """
        fastANI \\
            -q $query \\
            -r $reference \\
            -o ${prefix}.ani.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
        END_VERSIONS
        """
    }
}
