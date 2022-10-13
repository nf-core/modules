process HMMCOPY_MAPCOUNTER {
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_7':

    input:
    path bigwig

    output:
    path "*.map.wig"              , emit: wig
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '0.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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
