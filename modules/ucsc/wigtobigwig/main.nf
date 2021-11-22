def VERSION = '377' // No version information printed

process UCSC_WIGTOBIGWIG {
    tag '$wig'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ucsc-wigtobigwig=377" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-wigtobigwig:377--h0b8a92a_2' :
        'quay.io/biocontainers/ucsc-wigtobigwig:377--h0b8a92a_2' }"

    input:
    path wig
    path chromsizes

    output:
    path "*.bw"                   , emit: bw
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    wigToBigWig \\
        $args \\
        $wig \\
        $chromsizes \\
        ${wig.getSimpleName()}.bw

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        ucsc: \$(echo "$VERSION")
    END_VERSIONS
    """
}
