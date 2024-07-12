process PBTK_PBINDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbtk:3.1.1--h9ee0642_0':
        'biocontainers/pbtk:3.1.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.pbi"), emit: pbi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pbindex \\
        -j $task.cpus \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbindex: \$(pbindex --version | sed -n 's/pbindex \\(.*\\)/\\1/p')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.pbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbindex: \$(pbindex --version | sed -n 's/pbindex \\(.*\\)/\\1/p')
    END_VERSIONS
    """
}
