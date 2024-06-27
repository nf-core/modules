
process GMMDEMUX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gmm-demux:0.2.2.3--pyh7cba7a3_0':
        'biocontainers/gmm-demux:0.2.2.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(hto_matrix),val(hto_names)
    val type_report
    val summary_report
    val skip

    output:
    tuple val(meta), path("${meta.id}/barcodes.tsv.gz"                          ), emit: barcodes
    tuple val(meta), path("${meta.id}/*.mtx.gz"                                 ), emit: matrix
    tuple val(meta), path("${meta.id}/features.tsv.gz"                          ), emit: features
    tuple val(meta), path("${meta.id}/classification_report_${meta.id}"         ), emit: classification_report
    tuple val(meta), path("summary_report_${meta.id}.txt"                       ), emit: summary_report, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def type_report    = type_report ? "-f test/classification_report_${prefix}" : "-s test/classification_report_${prefix}"
    def skip           = skip ? "--skip $skip" : ""
    def summary_rep    = summary_report ? "-r summary_report_${prefix}.txt" : ""
    def VERSION        = '0.2.2.3' // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    """
    if [[ ${summary_report} == true ]]; then
        cat /dev/null >  summary_report_${prefix}.txt
        echo "summary report file created"
    fi 

    GMM-demux $args \\
        $type_report \\
        $summary_rep  \\
        $skip         \\
        $hto_matrix   \\
        $hto_names
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GMM-Demux: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.2.2.3'
    """
    mkdir test
    touch test/barcodes.tsv.gz
    touch test/features.tsv.gz
    touch test/matrix.mtx.gz
    touch test/classification_report_test

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GMM-Demux: $VERSION
    END_VERSIONS
    """
}
