process MINIBWA_INDEX {
    tag "$fasta"
    // NOTE minibwa builds an FM-index with libsais; peak memory scales with the reference size.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65579835ae1cda4cac4bec86b6210473b95f88d2ce8e4d7b6257a81bee1e8a27/data'
        : 'community.wave.seqera.io/library/minibwa:0.2--38d72f7ba96c74d9'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("minibwa"), emit: index
    tuple val("${task.process}"), val('minibwa'), eval('minibwa version | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_minibwa, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    """
    mkdir minibwa
    minibwa \\
        index \\
        $args \\
        -t ${task.cpus} \\
        $fasta \\
        minibwa/${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir minibwa
    touch minibwa/${prefix}.l2b
    touch minibwa/${prefix}.mbw
    """
}
