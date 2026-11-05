// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00cad0bb4c7cb85b0906f357e5197bcabb32569f9b48ff049f081b5bafd04aca/data'
        : 'community.wave.seqera.io/library/ncbi-amrfinderplus:4.2.7--01f7d2a8fad0d3aa'}"

    output:
    path "amrfinderdb.tar.gz", emit: db
    tuple val("${task.process}"), val('amrfinder'), eval('amrfinder --version'), emit: versions_amrfinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    amrfinder_update -d amrfinderdb
    tar czvf amrfinderdb.tar.gz -C amrfinderdb/\$(readlink amrfinderdb/latest) ./
    """

    stub:
    """
    touch amrfinderdb.tar
    gzip amrfinderdb.tar
    """
}
