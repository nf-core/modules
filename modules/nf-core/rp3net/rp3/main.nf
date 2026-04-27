process RP3NET_RP3 {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/25/2593eb0fc926d1d4c319576dc423fe7c7b9e4632638043560b927123e25b0e33/data' :
        'community.wave.seqera.io/library/pip_python_rp3net:a185caa894da5c7b' }"

    input:
    tuple val(meta), path(fasta), path(checkpoint)

    output:
    tuple val(meta), path("${prefix}.csv"), emit: scores
    tuple val("${task.process}"), val('rp3net'), eval("pip show RP3Net 2>&1 | awk '/^Version/{print \$2}'"), topic: versions, emit: versions_rp3net

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
