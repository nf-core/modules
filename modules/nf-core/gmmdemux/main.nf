
process GMMDEMUX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gmm-demux:0.2.2.3--pyh7cba7a3_0':
        'biocontainers/gmm-demux:0.2.2.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(cell_hashing_barcodes,stageAs: "hto_files/*"), path(cell_hashing_matrix,stageAs: "hto_files/*"),path(cell_hashing_features,stageAs: "hto_files/*"),val(hto_names)
    val type_report
    val skip

    output:
    tuple val(meta), path('test/barcodes.tsv.gz'), emit: barcodes
    tuple val(meta), path('test/*.mtx.gz'       ), emit: matrix
    tuple val(meta), path('test/features.tsv.gz'), emit: features
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def type_report    = type_report ? "-s simplfied_report_${prefix}" : "-f full_report_${prefix}"
    def skip           = skip ? "--skip $skip" : ""
    def VERSION        = '0.2.2.3' // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    """
    GMM-demux $args \\
        $type_report \\
        $skip \\
        hto_files \\
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


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GMM-Demux: $VERSION
    END_VERSIONS
    """
}
