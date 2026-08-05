process QCATCH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a7/a7d0112866550e3bcf97c40104596a3ca2ecbc26c13cf919fe76587554528281/data':
        'community.wave.seqera.io/library/pip_qcatch:03b88593a5cca75b' }"
    input:
    tuple val(meta), val(chemistry), path(quant_dir)

    output:
    tuple val(meta), path("*.html")                , emit: report
    tuple val(meta), path("*_filtered_quants.h5ad") , emit: filtered_h5ad
    tuple val(meta), path("*_metrics_summary.csv")  , emit: metrics_summary
    tuple val("${task.process}"), val('qcatch'), eval("qcatch --version | sed -e 's/qcatch //g'"), emit: versions_qcatch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chemistry_arg = chemistry ? "--chemistry ${chemistry}" : ''

    """
    qcatch \\
        --input ${quant_dir} \\
        --output ${prefix} \\
        ${chemistry_arg} \\
        --save_filtered_h5ad \\
        --export_summary_table \\
        ${args}

    mv ${prefix}/QCatch_report.html ${prefix}_qcatch_report.html
    mv ${prefix}/filtered_quants.h5ad ${prefix}_filtered_quants.h5ad
    mv ${prefix}/summary_table.csv ${prefix}_metrics_summary.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_qcatch_report.html
    touch ${prefix}_filtered_quants.h5ad
    touch ${prefix}_metrics_summary.csv
    """
}
