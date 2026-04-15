process WISECONDORX_GENDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/13af39819608398807612090d4b8af7dedb8db403967e71af22dbbeeb502ead1/data':
        'community.wave.seqera.io/library/wisecondorx:1.3.0--835c946afbce9082' }"

    input:
    tuple val(meta), path(npz)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), stdout, emit: gender
    tuple val("${task.process}"), val('wisecondorx'), eval("pip list |& sed -n 's/wisecondorx *//p'"), emit: versions_wisecondorx, topic: versions

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
