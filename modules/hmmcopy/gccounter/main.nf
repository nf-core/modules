process HMMCOPY_GCCOUNTER {
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_7' :
        'quay.io/biocontainers/hmmcopy:0.1.1--h2e03b76_7' }"

    input:
    path fasta

    output:
    path "*.gc.wig"    , emit: wig
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '0.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gcCounter \\
        $args \\
        ${fasta} > ${fasta.baseName}.gc.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmcopy: $VERSION
    END_VERSIONS
    """
}
