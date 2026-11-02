// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ANGSD_CONTAMINATION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a0185040abdace3529e7c8ae31a80bd86a47577e0d2078b2223dc1b1c362e656/data'
        : 'community.wave.seqera.io/library/angsd_htslib:714d1c4bd0020e05'}"

    input:
    tuple val(meta), path(icounts)
    tuple val(meta2), path(hapmap_file)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('angsd'), eval("angsd 2>&1 | sed '1!d;s/.*version: //;s/ .*//'"), emit: versions_angsd, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def seed_cmd = args.contains("-s ") ? '' : '-s 1'
    """
    contamination \
        ${args} \
        ${seed_cmd} \
        -a ${icounts} \
        -h ${hapmap_file} \
        -p ${task.cpus} \
        2>| >(tee ${prefix}.txt >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
