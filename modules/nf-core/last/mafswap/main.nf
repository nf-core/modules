process LAST_MAFSWAP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0d/0d27a2649f1291ff817dc8f73357ffac206424cd972d3855421e4258acc600f7/data'
        : 'community.wave.seqera.io/library/last:1611--e1193b3871fa0975'}"

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
    set -o pipefail
    maf-swap $args $maf | gzip --no-name > ${prefix}.swapped.maf.gz

    # maf-swap has no --version option but lastdb, part of the same package, has.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version 2>&1 | sed 's/lastdb //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub | gzip --no-name > ${prefix}.swapped.maf.gz

    # maf-swap has no --version option but lastdb, part of the same package, has.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version 2>&1 | sed 's/lastdb //')
    END_VERSIONS
    """
}
