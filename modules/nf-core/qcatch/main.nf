process QCATCH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/qcatch:0.2.8--3089b62e628f96d7':
        'community.wave.seqera.io/library/qcatch:0.2.8--454a9b478b62c36f' }"
    input:
    tuple val(meta), val(chemistry), path(quant_dir)

    output:
    tuple val(meta), path("*.html")                         , emit: report
    tuple val(meta), path("*_filtered_quants.h5ad") , emit: filtered_h5ad
    tuple val(meta), path("*_metrics_summary.csv")  , emit: metrics_summary
    path  "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export MPLCONFIGDIR=./tmp
    export NUMBA_CACHE_DIR=./tmp
    export NUMBA_DISABLE_CACHE=1
    export XDG_CACHE_HOME=./tmp
    export TMPDIR=./tmp

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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export MPLCONFIGDIR=./tmp
    export NUMBA_CACHE_DIR=./tmp
    export NUMBA_DISABLE_CACHE=1
    export XDG_CACHE_HOME=./tmp
    export TMPDIR=./tmp

    touch ${prefix}_qcatch_report.html
    touch ${prefix}_filtered_quants.h5ad
    touch ${prefix}_metrics_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qcatch: \$(qcatch --version | sed -e "s/qcatch //g")
    END_VERSIONS
    """
}
