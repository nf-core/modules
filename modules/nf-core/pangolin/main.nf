process PANGOLIN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pangolin=4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pangolin:4.2--pyhdfd78af_1' :
        'biocontainers/pangolin:4.2--pyhdfd78af_1' }"

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
