process METHYLSIEVE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92a979cffe85f24353179c5647b0e6151bf0c6d3bdbe8890101b6999107c4db5/data'
        : 'community.wave.seqera.io/library/methylsieve:0.1.0--dd751f2cc2a188b1'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("*.stats.tsv"), emit: stats
    tuple val("${task.process}"), val('methylsieve'), eval("methylsieve --version | sed 's/^methylsieve //'"), emit: versions_methylsieve, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    methylsieve \\
        --reference ${fasta} \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        --stats ${prefix}.stats.tsv \\
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
    touch ${prefix}.stats.tsv
    """
}
