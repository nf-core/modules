process DATAVZRD {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7e1ff6edeffa771f340c0c044fde514b61bc4691ede6b79728263bf20c7990f6/data'
        : 'community.wave.seqera.io/library/datavzrd:2.61.7--3dfa5861de05f6ac'}"

    input:
    tuple val(meta), path(config_file), path(table)

    output:
    tuple val(meta), path("${prefix}"), emit: report
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    datavzrd \\
        ${args} \\
        ${config_file} \\
        --output ${prefix} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datavzrd: \$(echo \$( datavzrd --version | sed -e 's/[^0-9.]//g' ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    mkdir ${prefix}/static
    mkdir ${prefix}/network
    mkdir ${prefix}/network/data
    mkdir ${prefix}/network/plots
    touch ./${prefix}/index.html
    touch ./${prefix}/static/bootstrap.min.css
    touch ./${prefix}/static/bootstrap-select.min.css
    touch ./${prefix}/static/bootstrap-table.min.css
    touch ./${prefix}/static/bootstrap-table-fixed-columns.min.css
    touch ./${prefix}/static/bundle.js
    touch ./${prefix}/static/datavzrd.css
    touch ./${prefix}/network/index_1.html
    touch ./${prefix}/network/config.js
    touch ./${prefix}/network/functions.js
    touch ./${prefix}/network/heatmap.js
    touch ./${prefix}/network/data/data_1.js
    touch ./${prefix}/network/plots/plot_0.js

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datavzrd: \$(echo \$( datavzrd --version | sed -e 's/[^0-9.]//g' ))
    END_VERSIONS
    """
}
