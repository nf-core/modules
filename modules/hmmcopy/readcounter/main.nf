def VERSION = '0.1.1' // Version information not provided by tool on CLI

process HMMCOPY_READCOUNTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_5' :
        'quay.io/biocontainers/hmmcopy:0.1.1--h2e03b76_5' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.wig"), emit: wig
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    readCounter \\
        $args \\
        ${bam} > ${prefix}.wig

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        hmmcopy: $VERSION
    END_VERSIONS
    """
}
