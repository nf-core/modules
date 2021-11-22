process LAST_SPLIT {
    tag "$meta.id"
    label 'process_high'

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
    """
    zcat < $maf | last-split $args | gzip --no-name > ${prefix}.maf.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(last-split --version 2>&1 | sed 's/last-split //')
    END_VERSIONS
    """
}
