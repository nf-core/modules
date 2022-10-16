process MERYL_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::meryl=1.3" : null)
    def container_image = "meryl:1.3--h87f3376_1"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.meryldb"), emit: meryl_db
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for READ in $reads; do
        meryl count \\
            threads=$task.cpus \\
            $args \\
            $reads \\
            output read.\${READ%.f*}.meryldb
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """
}
