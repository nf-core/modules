// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process WISECONDORX_GENDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/4103a257615e467895ebb0ff8aca0c17f80e852466feb36736236a4ca2597abf/data'
        : 'community.wave.seqera.io/library/wisecondorx:1.3.3--dd93194eb6a11e3b'}"

    input:
    tuple val(meta), path(npz)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), stdout, emit: gender
    tuple val("${task.process}"), val('wisecondorx'), eval("python -c \"import wisecondorx; print(wisecondorx.__version__)\""), emit: versions_wisecondorx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    WisecondorX gender \\
        ${npz} \\
        ${reference}
    """

    stub:
    """
    echo male
    """
}
