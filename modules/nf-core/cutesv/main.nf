process CUTESV {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4e499f094cacf5729232e96b56c3f13b770e2153e19d45e4ca3093969fafa22e/data' :
        'community.wave.seqera.io/library/cutesv:2.1.4--f0ec240b872b2a98' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val("cuteSV"), eval("cuteSV --version 2>&1 | sed 's/cuteSV //'"), topic: versions, emit: versions_cutesv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cuteSV \
        ${bam} \\
        ${fasta} \\
        ${prefix}.vcf \\
        . \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
