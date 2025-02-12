process LAST_POSTMASK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db0b5de918238f07ec1ca668be942397da85e26aa582f8927ac37c70896303cf/data'
        : 'community.wave.seqera.io/library/last:1608--f41c047f7dc37e30'}"

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
    set -o pipefail
    last-postmask $args $maf | gzip --no-name > ${prefix}.maf.gz

    # last-postmask does not have a --version option
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastal --version 2>&1 | sed 's/lastal //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo stub | gzip --no-name > ${prefix}.maf.gz

    # last-postmask does not have a --version option
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastal --version 2>&1 | sed 's/lastal //')
    END_VERSIONS
    """
}
