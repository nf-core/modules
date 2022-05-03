def VERSION = '0.1.1'

process HMMCOPY_GENERATEMAP {
    tag '$bam'
    label 'process_long'

    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_7':
        'quay.io/biocontainers/hmmcopy:0.1.1--h2e03b76_7' }"

    input:
    path fasta

    output:
    path "*.map.bw"              , emit: bigwig
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

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
