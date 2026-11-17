// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FFQ {
    tag "${ids.size() == 1 ? ids[0] : "${ids[0]+'..'+ids[-1]}"}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f651420b3d10407351f47545bb8d274e3be9d98c709a348696de84101e36a95c/data'
        : 'community.wave.seqera.io/library/ffq_python:1b4df44896918bd7'}"

    input:
    val ids

    output:
    path "*.json"                                                                                              , emit: json
    tuple val("${task.process}"), val('ffq'), eval("ffq --help 2>&1 | sed 's/^.*ffq //; s/: A command.*\$//'") , emit: versions_ffq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def id_list = ids.sort()
    def name = id_list.size() == 1 ? ids[0] : 'metadata'
    def prefix = task.ext.prefix ?: "${name}"
    """
    ffq \\
        ${id_list.join(' ')} \\
        $args \\
        > ${prefix}.json
    """

    stub:
    def id_list = ids.sort()
    def name = id_list.size() == 1 ? ids[0] : 'metadata'
    def prefix = task.ext.prefix ?: "${name}"
    """
    touch ${prefix}.json
    """
}
