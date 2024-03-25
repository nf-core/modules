process GTFSORT {
    tag "$gtf"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtfsort:0.2.2--h4ac6f70_0':
        'biocontainers/gtfsort:0.2.2--h4ac6f70_0' }"

    input:

    path gtf

    output:
    path "*.sorted.gtf", emit: gtf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"

    """
    gtfsort \\
        -i $gtf \\
        -o ${prefix}.sorted.gtf \\
        -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtfsort: \$(gtfsort --version |& sed 's/gtfsort //')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"

    """
    touch ${prefix}.sorted.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtfsort: \$(gtfsort --version |& sed 's/gtfsort //')
    END_VERSIONS
    """
}
