process QCATCH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/6233ebf32d1ce85d41398121c4b8c1449779fa13b1529e702e067bfca6f835b7/data':
        'community.wave.seqera.io/library/pip_qcatch:05d150ed59c3ae84' }"
    input:
    tuple val(meta), val(chemistry), path(quant_dir)

    output:
    tuple val(meta), path("*.html")                , emit: report
    tuple val(meta), path("*_filtered_quants.h5ad") , emit: filtered_h5ad
    tuple val(meta), path("*_metrics_summary.csv")  , emit: metrics_summary
    path  "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    qcatch \\
        --input ${quant_dir} \\
        --output ${prefix} \\
        --chemistry ${chemistry} \\
        --save_filtered_h5ad \\
        --export_summary_table \\
        ${args}

    mv ${prefix}/QCatch_report.html ${prefix}_qcatch_report.html
    mv ${prefix}/filtered_quants.h5ad ${prefix}_filtered_quants.h5ad
    mv ${prefix}/summary_table.csv ${prefix}_metrics_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qcatch: \$(qcatch --version | sed -e "s/qcatch //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_qcatch_report.html
    touch ${prefix}_filtered_quants.h5ad
    touch ${prefix}_metrics_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qcatch: \$(qcatch --version | sed -e "s/qcatch //g")
    END_VERSIONS
    """
}
