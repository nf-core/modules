process PLOTSR {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/plotsr:1.1.1--pyh7cba7a3_0'
        : 'biocontainers/plotsr:1.1.1--pyh7cba7a3_0'}"

    input:
    tuple val(meta), path(syri)
    tuple val(meta2), path(fastas)
    tuple val(meta3), path(genomes)
    tuple val(meta4), path(bedpe)
    tuple val(meta5), path(markers)
    tuple val(meta6), path(tracks)
    tuple val(meta7), path(chrord)
    tuple val(meta8), path(chrname)

    output:
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def syri_arg = syri instanceof List ? syri.collect { syri_file -> "--sr ${syri_file}" }.join(' ') : "--sr ${syri}"
    def bedpe_arg = bedpe ? "--bedpe ${bedpe}" : ''
    def markers_arg = markers ? "--markers ${markers}" : ''
    def tracks_arg = tracks ? "--tracks ${tracks}" : ''
    def chrord_arg = chrord ? "--chrord ${chrord}" : ''
    def chrname_arg = chrname ? "--chrname ${chrname}" : ''
    """
    plotsr \\
        ${syri_arg} \\
        --genomes ${genomes} \\
        ${bedpe_arg} \\
        ${markers_arg} \\
        ${tracks_arg} \\
        ${chrord_arg} \\
        ${chrname_arg} \\
        ${args} \\
        -o ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotsr: \$(plotsr --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotsr: \$(plotsr --version)
    END_VERSIONS
    """
}
