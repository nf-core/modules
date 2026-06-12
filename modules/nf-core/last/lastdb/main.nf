process LAST_LASTDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2e/2eb57450207840a7fba7f60b65239a86679bfcaa79fb5fba652dd41af2b3e1d9/data'
        : 'community.wave.seqera.io/library/last:1651--c83f04148c23181f'}"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("lastdb"), emit: index
    // last-dotplot has no --version option so let's use lastal from the same suite
    tuple val("${task.process}"), val('last'), eval("lastal --version | sed 's/lastal //'"), emit: versions_last, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir lastdb
    lastdb \\
        $args \\
        -P $task.cpus \\
        lastdb/${prefix} \\
        $fastx
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir lastdb
    touch lastdb/${prefix}.bck
    touch lastdb/${prefix}.des
    touch lastdb/${prefix}.prj
    touch lastdb/${prefix}.sds
    touch lastdb/${prefix}.ssp
    touch lastdb/${prefix}.suf
    touch lastdb/${prefix}.tis

    """
}
