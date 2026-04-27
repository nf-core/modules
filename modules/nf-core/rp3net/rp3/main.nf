process RP3NET_RP3 {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/1653ee5e71a6b4ff66690a7b18238aae9b5b78dd38b9b4aed58057e06693d717/data' :
        'community.wave.seqera.io/library/pip_python_rp3net_peft:99495b7509bc1482' }"

    input:
    tuple val(meta), path(fasta), path(checkpoint)

    output:
    tuple val(meta), path("${prefix}.csv"), emit: scores
    tuple val("${task.process}"), val('rp3net'), eval("python -c \"import importlib.metadata; print(importlib.metadata.version('RP3Net'))\""), topic: versions, emit: versions_rp3net

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    rp3 \\
        -p ${checkpoint} \\
        -f ${fasta} \\
        -o ${prefix}.csv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
