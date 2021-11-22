process LAST_MAFSWAP {
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
    """
    maf-swap $args $maf | gzip --no-name > ${prefix}.swapped.maf.gz

    # maf-swap has no --version option but lastdb, part of the same package, has.
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(lastdb --version 2>&1 | sed 's/lastdb //')
    END_VERSIONS
    """
}
