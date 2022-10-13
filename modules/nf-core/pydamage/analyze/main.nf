process PYDAMAGE_ANALYZE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pydamage=0.70" : null)
    def container_image = "/pydamage:0.70--pyhdfd78af_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("pydamage_results/pydamage_results.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pydamage \\
        analyze \\
        $args \\
        -p $task.cpus \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pydamage: \$(echo \$(pydamage --version 2>&1) | sed -e 's/pydamage, version //g')
    END_VERSIONS
    """
}
