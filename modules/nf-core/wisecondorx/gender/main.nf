process WISECONDORX_GENDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/10a1cb134d692ea5d26cb7b4a5a2e83fe28eacac5b637a8eb6ca4d9532602222/data':
        'community.wave.seqera.io/library/wisecondorx:1.3.2--7ece9fb804446823' }"

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
