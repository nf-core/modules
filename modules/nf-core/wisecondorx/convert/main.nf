process WISECONDORX_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/10a1cb134d692ea5d26cb7b4a5a2e83fe28eacac5b637a8eb6ca4d9532602222/data':
        'community.wave.seqera.io/library/wisecondorx:1.3.2--7ece9fb804446823' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.npz"), emit: npz
    tuple val("${task.process}"), val('wisecondorx'), eval("python -c \"import wisecondorx; print(wisecondorx.__version__)\""), emit: versions_wisecondorx, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""

    """
    WisecondorX convert \\
        ${bam} \\
        ${prefix} \\
        ${reference} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.npz
    """
}
