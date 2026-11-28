// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SINGLEM_DBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c64fb959b73c0d5cefcd3fff9c0d8f388f635ab1e953cb658adf267097eb77df/data'
        : 'community.wave.seqera.io/library/singlem:0.21.4--12242454960ab14c'}"

    output:
    path("*.smpkg.zb")                                                                                                     , emit: singlem_database
    tuple val("${task.process}"), val('singlem'), eval('singlem --version'), topic: versions                                , emit: versions_singlem

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    singlem \\
        data \\
        --output-directory . \\
        ${args}
    """

    stub:
    """
    touch S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb
    """
}
