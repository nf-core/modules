process DATAVZRD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/datavzrd:2.23.2':
        'biocontainers/datavzrd:2.23.2' }"

    input:
    tuple val(meta), file (config_file)

    output:
    tuple val(meta), path ("output"), emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    mkdir output
    datavzrd \\
        ${config_file} \\
        --output output \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datavzrd: \$(echo \$( datavzrd version | sed -e "s/datavzrd v//g" ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
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
        datavzrd: \$(echo \$( datavzrd version | sed -e "s/datavzrd v//g" ))
    END_VERSIONS
    """
}
