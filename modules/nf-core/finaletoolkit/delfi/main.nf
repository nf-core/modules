process FINALETOOLKIT_DELFI {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34f01d128ed135aedc33ddb62fced3911bef6d1a909694291b7184bf83719402/data'
        : 'community.wave.seqera.io/library/finaletoolkit:0.11.1--8fe5ba6ec9e2ec95'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(genome_2bit)
    tuple val(meta3), path(chromosome_sizes)
    tuple val(meta4), path(bins)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('finaletoolkit'), eval("finaletoolkit --version | sed 's/FinaleToolkit //g'"), topic: versions, emit: versions_finaletoolkit


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_delfi"
    """
    finaletoolkit \\
        delfi \\
        -w ${task.cpus} \\
        ${bam} \\
        ${chromosome_sizes} \\
        ${genome_2bit} \\
        ${bins} \\
        ${args} \\
        -o "${prefix}.bed"
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_delfi"
    """
    echo ${args}

    touch ${prefix}.bed
    """
}
