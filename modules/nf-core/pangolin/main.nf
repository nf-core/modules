process PANGOLIN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::pangolin=4.1.1' : null)
    def container_image = "/pangolin:4.1.1--pyhdfd78af_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.csv'), emit: report
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pangolin \\
        $fasta\\
        --outfile ${prefix}.pangolin.csv \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangolin: \$(pangolin --version | sed "s/pangolin //g")
    END_VERSIONS
    """
}
