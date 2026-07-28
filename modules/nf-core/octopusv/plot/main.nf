process OCTOPUSV_PLOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/octopusv:0.3.3--pyhdfd78af_0'
        : 'quay.io/biocontainers/octopusv:0.3.3--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path("${prefix}_sv_sizes.png"), emit: sv_sizes_png
    tuple val(meta), path("${prefix}_sv_sizes.svg"), emit: sv_sizes_svg
    tuple val(meta), path("${prefix}_sv_types.png"), emit: sv_types_png
    tuple val(meta), path("${prefix}_sv_types.svg"), emit: sv_types_svg
    tuple val(meta), path("${prefix}_chromosome_distribution.png"), emit: chromosome_distribution_png
    tuple val(meta), path("${prefix}_chromosome_distribution.svg"), emit: chromosome_distribution_svg
    tuple val("${task.process}"), val('octopusv'), eval("python -c \"import importlib.metadata as m; print(m.version('octopusv'))\""), emit: versions_octopusv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    octopusv plot \\
        ${txt} \\
        --output-prefix ${prefix} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}_sv_sizes.svg
    touch ${prefix}_sv_sizes.png
    touch ${prefix}_sv_types.png
    touch ${prefix}_sv_types.svg
    touch ${prefix}_chromosome_distribution.png
    touch ${prefix}_chromosome_distribution.svg
    """
}
