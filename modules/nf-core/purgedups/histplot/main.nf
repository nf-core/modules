process PURGEDUPS_HISTPLOT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::purge_dups=1.2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--py39h7132678_1':
        'biocontainers/purge_dups:1.2.6--py39h7132678_1' }"

    input:
    tuple val(meta), path(statfile), path(cutoff)

    output:
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hist_plot.py \\
        -c $cutoff \\
        $args \\
        $statfile \\
        ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hist_plot : \$( hist_plot.py -v | sed 's/hist_plot //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hist_plot : \$( hist_plot.py -v | sed 's/hist_plot //' )
    END_VERSIONS
    """
}
