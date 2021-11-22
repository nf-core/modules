def VERSION = '0.1.1'

process HMMCOPY_READCOUNTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::hmmcopy=0.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_5"
    } else {
        container "quay.io/biocontainers/hmmcopy:0.1.1--h2e03b76_5"
    }

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("*.wig"), emit: wig
        path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
        def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    readCounter \\
        $args \\
        ${bam} > ${prefix}.wig

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
