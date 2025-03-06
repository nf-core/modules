process BUSCO_GENERATEPLOT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.7.1--pyhdfd78af_0':
        'biocontainers/busco:5.7.1--pyhdfd78af_0' }"

    input:
    path short_summary_txt, stageAs: 'busco/*'

    output:
    path '*.png'        , emit: png
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: 'busco_figure'
    """
    generate_plot.py \\
        $args \\
        -wd busco

    mv ./busco/busco_figure.png ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: 'busco_figure'
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
