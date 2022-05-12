def VERSION = '0.1.1' // Version information not provided by tool on CLI

process HMMCOPY_MAPCOUNTER {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_7':
        'quay.io/biocontainers/hmmcopy:0.1.1--h2e03b76_7' }"

    input:
    path bigwig

    output:
    path "*.map.wig"              , emit: wig
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mapCounter \\
        $args \\
        $bigwig > ${bigwig.baseName}.map.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmcopy: \$(echo $VERSION)
    END_VERSIONS
    """
}
