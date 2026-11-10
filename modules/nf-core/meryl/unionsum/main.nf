// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MERYL_UNIONSUM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7df6a8db5e64dc18ecc1aa31b12308fc9d5e07ddefcc64f5ffffa1716c029da8/data'
        : 'community.wave.seqera.io/library/meryl:1.4.1--2481ef1dd13fafa3'}"

    input:
    tuple val(meta), path(meryl_dbs)
    val kvalue

    output:
    tuple val(meta), path("*.unionsum.meryl"), emit: meryl_db
    // meryl --version reports a git revision, not the release; keep in sync with environment.yml
    tuple val("${task.process}"), val('meryl'), val('1.4.1'), emit: versions_meryl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    meryl union-sum \\
        k=$kvalue \\
        threads=$task.cpus \\
        memory=${task.memory.toGiga()} \\
        $args \\
        output ${prefix}.unionsum.meryl \\
        $meryl_dbs
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unionsum.meryl
    """
}
