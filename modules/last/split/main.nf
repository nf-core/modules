process LAST_SPLIT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/last:1250--h2e03b76_0' :
        'quay.io/biocontainers/last:1250--h2e03b76_0' }"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat < $maf | last-split $args | gzip --no-name > ${prefix}.maf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(last-split --version 2>&1 | sed 's/last-split //')
    END_VERSIONS
    """
}
