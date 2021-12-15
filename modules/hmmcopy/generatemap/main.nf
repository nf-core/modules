process HMMCOPY_GENERATEMAP {
    tag '$bam'
    label 'process_long'

    def VERSION = '0.1.1'

    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_7':
        'quay.io/biocontainers/hmmcopy:0.1.1--h2e03b76_7' }"

    input:
    path fasta

    output:
    path "*.map.bw"              , emit: bigwig
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    # build required indexes
    perl /usr/local/bin/generateMap.pl -b \\
        $args \\
        $fasta

    # run
    perl /usr/local/bin/generateMap.pl \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmcopy: \$(echo $VERSION)
    END_VERSIONS
    """
}
