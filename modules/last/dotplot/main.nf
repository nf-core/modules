process LAST_DOTPLOT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/last:1250--h2e03b76_0' :
        'quay.io/biocontainers/last:1250--h2e03b76_0' }"

    input:
    tuple val(meta), path(maf)
    val(format)

    output:
    tuple val(meta), path("*.gif"), optional:true, emit: gif
    tuple val(meta), path("*.png"), optional:true, emit: png
    path "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    last-dotplot \\
        $args \\
        $maf \\
        $prefix.$format

    # last-dotplot has no --version option so let's use lastal from the same suite
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(lastal --version | sed 's/lastal //')
    END_VERSIONS
    """
}
