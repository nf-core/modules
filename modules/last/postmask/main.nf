process LAST_POSTMASK {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/last:1250--h2e03b76_0' :
        'quay.io/biocontainers/last:1250--h2e03b76_0' }"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use the suffix option to disambiguate"
    """
    last-postmask $args $maf | gzip --no-name > ${prefix}.maf.gz

    # last-postmask does not have a --version option
    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        last: \$(lastal --version 2>&1 | sed 's/lastal //')
    END_VERSIONS
    """
}
