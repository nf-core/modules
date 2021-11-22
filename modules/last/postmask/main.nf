process LAST_POSTMASK {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/last:1250--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/last:1250--h2e03b76_0"
    }

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    path "versions.yml"              , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use the suffix option to disambiguate"
    """
    last-postmask $args $maf | gzip --no-name > ${prefix}.maf.gz

    # last-postmask does not have a --version option
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(lastal --version 2>&1 | sed 's/lastal //')
    END_VERSIONS
    """
}
