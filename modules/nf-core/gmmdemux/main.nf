
process GMMDEMUX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gmm-demux:0.2.2.3--pyh7cba7a3_0':
        'biocontainers/gmm-demux:0.2.2.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(cell_hashing_barcodes,stageAs: "hto_files/*"), path(cell_hashing_matrix,stageAs: "hto_files/*"),path(cell_hashing_features,stageAs: "hto_files/*"),val(hto_names)
    val csv
    val simplified_report
    val examine
    val ambigous
    val extract
    val skip

    output:
    tuple val(meta), path('test/barcodes.tsv.gz'), emit: barcodes
    tuple val(meta), path('test/*.mtx.gz'       ), emit: matrix
    tuple val(meta), path('test/features.tsv.gz'), emit: features
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //since this tool has many optional inputs that can be passed in, we need to check if they are null or not
    // in order to produce or not certain reports or outputs
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def simplified_rep = simplified_report ? "--simplified $simplified_report" : ""
    def examine = examine ? "--examine $examine" : ""
    def extract = extract ? "--extract $extract" : ""
    def skip = skip ? "--skip $skip" : ""
    def type = csv ? "-c" : ""
    def ambigous = examine ? "--ambigous $ambigous " : ""
    def VERSION = '0.2.2.3' // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    """
    GMM-demux $type hto_files $hto_names $simplified_rep $examine $ambigous $extract $skip $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GMM-Demux: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
