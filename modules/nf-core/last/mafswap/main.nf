process LAST_MAFSWAP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2e/2eb57450207840a7fba7f60b65239a86679bfcaa79fb5fba652dd41af2b3e1d9/data'
        : 'community.wave.seqera.io/library/last:1651--c83f04148c23181f'}"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    // last-dotplot has no --version option so let's use lastal from the same suite
    tuple val("${task.process}"), val('last'), eval("lastal --version | sed 's/lastal //'"), emit: versions_last, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail
    maf-swap $args $maf | gzip --no-name > ${prefix}.swapped.maf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.swapped.maf.gz
    """
}
