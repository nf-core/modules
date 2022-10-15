process HMMCOPY_GENERATEMAP {
    tag "$bam"
    label 'process_long'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
    def container_image = "hmmcopy:0.1.1--h2e03b76_7"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    path fasta

    output:
    path "*.map.bw"              , emit: bigwig
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '0.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    # build required indexes
    generateMap.pl -b \\
        $args \\
        $fasta

    # run
    generateMap.pl \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmcopy: \$(echo $VERSION)
    END_VERSIONS
    """
}
