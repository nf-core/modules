process FOLDCOMP_COMPRESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldcomp:0.0.7--h43eeafb_0':
        'biocontainers/foldcomp:0.0.7--h43eeafb_0' }"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("*fcz"), emit: fcz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = pdb.isDirectory() ? "_fcz" : ".fcz"
    """
    foldcomp \\
        compress \\
        -t ${task.cpus} \\
        ${pdb} \\
        ${prefix}${ext} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldcomp: \$(foldcomp --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fcz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldcomp: \$(foldcomp --version)
    END_VERSIONS
    """
}
