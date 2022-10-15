process KRAKENTOOLS_COMBINEKREPORTS {
    label 'process_single'

    conda (params.enable_conda ? "bioconda::krakentools=1.2" : null)
    def container_image = "krakentools:1.2--pyh5e36f6f_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(kreports)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    combine_kreports.py \\
        -r ${kreports} \\
        -o ${prefix}.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combine_kreports.py: ${VERSION}
    END_VERSIONS
    """
}
