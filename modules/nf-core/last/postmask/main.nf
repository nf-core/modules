process LAST_POSTMASK {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    def container_image = "/last:1250--h2e03b76_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }
n

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
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    last-postmask $args $maf | gzip --no-name > ${prefix}.maf.gz

    # last-postmask does not have a --version option
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastal --version 2>&1 | sed 's/lastal //')
    END_VERSIONS
    """
}
