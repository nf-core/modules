process BUSCO_GENERATEPLOT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c607f319867d96a38c8502f751458aa78bbd18fe4c7c4fa6b9d8350e6ba11ebe/data'
        : 'community.wave.seqera.io/library/busco_sepp:f2dbc18a2f7a5b64'}"

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
        busco: \$( busco --version 2> /dev/null | sed 's/BUSCO //g' )
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: 'busco_figure'
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2> /dev/null | sed 's/BUSCO //g' )
    END_VERSIONS
    """
}
