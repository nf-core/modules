// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BEDOPS_GTF2BED {
    tag "$gtf"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46dc39b50ab3b00b757de451db624d1fdec87ee14a0cd5234598b71f200d92e0/data'
        : 'community.wave.seqera.io/library/bedops:2.4.42--1f31d470625cb9a3' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path('*.bed'), emit: bed
    tuple val("${task.process}"), val("bedops"), eval('bedops --version | sed -n "s/.*version: *\\([^ ]*\\).*/\\1/p"'), emit: versions_bedops, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"

    """
    cat \\
    $gtf \\
    | gtf2bed \\
    $args \\
    --attribute-key=exon_id \\
    > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    touch ${prefix}.bed
    """

}
