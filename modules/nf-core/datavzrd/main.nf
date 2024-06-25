process DATAVZRD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/datavzrd:2.36.12--bb93c8c988b7a9af':
        'community.wave.seqera.io/library/datavzrd:2.36.12--593eb75e566b7f2a' }"

    input:
    tuple val(meta), file (config_file), file(table)

    output:
    tuple val(meta), path ("output"), emit: report
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir output
    datavzrd \\
        ${args} \\
        ${config_file} \\
        --output output \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datavzrd: \$(echo \$( datavzrd --version | sed -e 's/[^0-9.]//g' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir output
    mkdir output/static
    mkdir output/network
    mkdir output/network/data
    mkdir output/network/plots
    touch ./output/index.html
    touch ./output/static/bootstrap.min.css
    touch ./output/static/bootstrap-select.min.css
    touch ./output/static/bootstrap-table.min.css
    touch ./output/static/bootstrap-table-fixed-columns.min.css
    touch ./output/static/bundle.js
    touch ./output/static/datavzrd.css
    touch ./output/network/index_1.html
    touch ./output/network/config.js
    touch ./output/network/functions.js
    touch ./output/network/heatmap.js
    touch ./output/network/data/data_1.js
    touch ./output/network/plots/plot_0.js

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datavzrd: \$(echo \$( datavzrd --version | sed -e 's/[^0-9.]//g' ))
    END_VERSIONS
    """
}
