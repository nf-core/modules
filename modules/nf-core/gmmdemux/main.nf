
process GMMDEMUX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gmm-demux:0.2.2.3--pyh7cba7a3_0':
        'biocontainers/gmm-demux:0.2.2.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(hto_matrix), val(hto_names)
    val type_report
    val summary_report
    path skip
    path examine

    output:
    tuple val(meta), path("barcodes.tsv.gz"     ), emit: barcodes
    tuple val(meta), path("matrix.mtx.gz"       ), emit: matrix
    tuple val(meta), path("features.tsv.gz"     ), emit: features
    tuple val(meta), path("GMM_*.csv"           ), emit: classification_report
    tuple val(meta), path("GMM_*.config"        ), emit: config_report
    tuple val(meta), path("summary_report_*.txt"), emit: summary_report, optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def skip           = skip ? "--skip $skip" : ""
    def examine_cells  = examine ? "--examine $examine" : ""
    def VERSION        = '0.2.2.3' // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    def type_report    = type_report ? "-f ." : "-s ."
    def summary_rep    = summary_report ? "-r ${prefix}_summary_report.txt" : ""
    """
    if [[ ${summary_report} == true ]]; then
        cat /dev/null > ${prefix}_summary_report.txt
    fi

    GMM-demux $args \\
        $type_report \\
        $summary_rep \\
        $skip \\
        $examine_cells \\
        $hto_matrix \\
        $hto_names \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GMM-Demux: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '0.2.2.3'
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > barcodes.tsv.gz
    echo "" | gzip > features.tsv.gz
    echo "" | gzip > matrix.mtx.gz
    touch GMM_full.config
    touch GMM_full.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GMM-Demux: $VERSION
    END_VERSIONS
    """
}
