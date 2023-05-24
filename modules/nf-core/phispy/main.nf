process PHISPY {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::phispy=4.2.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phispy:4.2.21--py310h30d9df9_1':
        'biocontainers/phispy:4.2.21--py310h30d9df9_1' }"

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("${prefix}/")                        , emit: results
    tuple val(meta), path("${prefix}/prophage_coordinates.tsv"), emit: coordinates
    tuple val(meta), path("${prefix}/*.gb*")                   , emit: gbk
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    PhiSpy.py \\
        $args \\
        --threads $task.cpus \\
        -o $prefix \\
        $gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PhiSpy: \$(echo \$(PhiSpy.py --version 2>&1))
    END_VERSIONS
    """
}
