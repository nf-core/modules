process PERBASE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perbase:0.10.2--h15397dd_0':
        'biocontainers/perbase:0.10.2--h15397dd_0' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    perbase \\
        base-depth \\
        $bam \\
        $args \\
        --threads $task.cpus \\
        --bgzip \\
        --output ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perbase: \$(perbase --version |& sed '1!d ; s/perbase //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perbase: \$(perbase --version |& sed '1!d ; s/perbase //')
    END_VERSIONS
    """
}
