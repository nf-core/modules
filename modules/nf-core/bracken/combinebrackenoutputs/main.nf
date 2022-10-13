process BRACKEN_COMBINEBRACKENOUTPUTS {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bracken=2.7" : null)
    def container_image = "/bracken:2.7--py39hc16433a_0"
                                                                container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    path input

    output:
    path "*.txt"       , emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "bracken_combined.txt"
    // WARN: Version information not provided by tool on CLI.
    // Please update version string below when bumping container versions.
    def VERSION = '2.7'
    """
    combine_bracken_outputs.py \\
        $args \\
        --files ${input} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combine_bracken_output: ${VERSION}
    END_VERSIONS
    """
}
