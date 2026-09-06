process DUPBLASTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ac/ac756f5169c3a6e6e0ab5884bbdc49171a6d0f669066b7494df75185ad4b0716/data'
        : 'community.wave.seqera.io/library/dupblaster:0.1.1--e033684cd2346e47'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("*.tsv"), emit: stats
    tuple val("${task.process}"), val('dupblaster'), eval("dupblaster --version | sed 's/^dupblaster //'"), emit: versions_dupblaster, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    dupblaster \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --stats ${prefix}.dupblaster.tsv \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    touch ${prefix}.dupblaster.tsv
    """
}
