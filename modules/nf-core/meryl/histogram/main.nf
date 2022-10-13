process MERYL_HISTOGRAM {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::meryl=1.3" : null)
        def container_image = "/meryl:1.3--h87f3376_1"
                                                  container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(meryl_db)

    output:
    tuple val(meta), path("*.hist"), emit: hist
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    meryl histogram \\
        threads=$task.cpus \\
        $args \\
        $meryl_db > ${prefix}.hist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """
}
