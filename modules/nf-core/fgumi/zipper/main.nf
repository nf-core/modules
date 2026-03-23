process FGUMI_ZIPPER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/954170443a820787c9e02ef2135ebb8ec29c6b03633b0d61b5fafa98c59a1cce/data'
        : 'community.wave.seqera.io/library/fgumi_r-base_r-ggplot2_r-scales:09c99070b82c1c28'}"

    input:
    tuple val(meta), path(mapped, stageAs: "mapped/*"), path(unmapped, stageAs: "unmapped/*")
    tuple val(meta2), path(fasta), path(fasta_index), path(fasta_dict)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval("fgumi --version | sed -n 's/^fgumi //p'"), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fgumi zipper \\
        ${args} \\
        --input ${mapped} \\
        --unmapped ${unmapped} \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        --output ${prefix}.bam
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.bam
    """
}
