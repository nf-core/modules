process KRAKENTOOLS_COMBINEKREPORTS {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::krakentools=1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
        'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0' }"

    input:
    path kreports

    output:
    path "*.txt", emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "kreports_combined.txt"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    combine_kreports.py \\
        -r ${kreports} \\
        -o ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combine_kreports.py: ${VERSION}
    END_VERSIONS
    """
}
